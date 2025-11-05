library(tidyverse)
library(ggpubr)
library(ggprism)
library(Seurat)
conflicted::conflict_prefer_all("dplyr", quiet = TRUE)

config_dir <- "./config/"
data_dir <- "../data_objects/"
figures_dir <- "../figures/"

# Parallell Processing
if(parallel::detectCores(logical=FALSE) > 3) {
  library(doParallel)
  
  num_cores <- parallel::detectCores(logical=FALSE)
  cl <- makePSOCKcluster(num_cores - 2)
  registerDoParallel(cl)
}

source("OmicsToolbox/marker_genes.R")
source("OmicsToolbox/seurat_norm_integrate.R")

##########################################

## Run the 10x Objects through the Seurat Pipeline ##

# configs <- c("sept_cornea_viral.yaml")
# 
# for (config in configs) {
#   loaded_config <- yaml::read_yaml(paste0(config_dir, config))
# 
#   seurat_norm_integrate(loaded_config, seurat_verbose = TRUE)
# }


##########################################

## Annotate the Cornea object ##

# Read in Seurat Object
sept_cornea <- SeuratDisk::LoadH5Seurat(paste0(data_dir, "seurat/sept_cornea_viral.h5Seurat"))

## Factorizations
sept_cornea$condition <- factor(sept_cornea$condition, levels = c("kos", "re"))
sept_cornea$id <- factor(sept_cornea$id)
sept_cornea$date <- factor(sept_cornea$date)
sept_cornea$location <- factor(sept_cornea$location)
sept_cornea$celltype <- factor(sept_cornea$celltype)


##########################################
VlnPlot(sept_cornea, assay = "SCT", 
        features = marker_genes_L2, 
        group.by = "seurat_clusters",
        #cols = scCustomize::DiscretePalette_scCustomize(length(marker_genes_L2), palette = "glasbey", shuffle = FALSE),
        stack = TRUE, flip = TRUE) + NoLegend() + 
  ggtitle("Sept Cornea")

sept_cornea_L2_avgs <- AverageExpression(sept_cornea, assays = "SCT", layer = "scale.data", group.by = "seurat_clusters",
                                         features = c("Ptprc", "Epcam", "Itgam", "Cd74", "Cd68", "Cd3g", "Pax5", "Klrb1c", "Vwf", "Ncam1", "Pdgfra"))[[1]] %>% 
  as.matrix() %>% t() %>% 
  as_tibble(rownames = "cluster") %>% 
  mutate(designation = case_when(Epcam > 0.4 ~ "Epi", 
                                 Cd3g > 1 ~ "T", 
                                 (Itgam > 0) | ((Cd74 > 0) & (Pax5 < 1))  ~ "Mye", 
                                 Klrb1c > 0 ~ "NK", 
                                 Pax5 > 1 ~ "B", 
                                 (Ptprc > 0) | (Cd74> 0) ~ "Mye",
                                 Vwf > 1 ~ "Endo", 
                                 Ncam1 > 1 & Pdgfra > 1 ~ "Ncam1+Pdgfra+",
                                 .default = "other"))

print(sept_cornea_L2_avgs, n = Inf)

sept_cornea$cell_L2 <- factor(sept_cornea_L2_avgs$designation[as.numeric(sept_cornea$seurat_clusters)], 
                              levels = c("Mye", "NK", "B", "T", "Epi", "Endo", "Ncam1+Pdgfra+", "other"))

sept_cornea$cell_L1 <- 
  factor(case_when(
    "Mye" == sept_cornea$cell_L2 ~ "immune", 
    "T" == sept_cornea$cell_L2 ~ "immune", 
    "B" == sept_cornea$cell_L2 ~ "immune", 
    "NK" == sept_cornea$cell_L2 ~ "immune", 
    "other_imm" == sept_cornea$cell_L2 ~ "immune", 
    "Epi" == sept_cornea$cell_L2 ~ "epithelial", 
    .default = "other"
  ), levels = c("immune", "epithelial", "other"))

## Cell L4 annotation
sept_cornea_cell_L4_mappings <- 
  sept_cornea@meta.data %>% 
  select(cell_L2, seurat_clusters) %>% 
  arrange(cell_L2, seurat_clusters) %>% 
  group_by(cell_L2) %>% 
  distinct() %>% 
  mutate(idx = row_number() - 1, 
         cell_L4 = paste0(cell_L2, "_", idx)) %>% 
  select(-idx) %>% 
  right_join(sept_cornea@meta.data %>% select(cell_L2, seurat_clusters) %>% rownames_to_column("barcode"), 
             join_by("cell_L2" == "cell_L2", "seurat_clusters" == "seurat_clusters")) %>% 
  column_to_rownames("barcode") %>% 
  select(cell_L4)

sept_cornea <- AddMetaData(sept_cornea, sept_cornea_cell_L4_mappings)
sept_cornea$cell_L4 <- factor(sept_cornea$cell_L4, 
                              levels =  str_sort(unique(sept_cornea$cell_L4), numeric = TRUE))


## Write the new annotated version and delete our unannotated version
sept_cornea %>% SeuratDisk::SaveH5Seurat(paste0(data_dir, "seurat/sept_cornea_viral_annot.h5Seurat"), overwrite = TRUE)

file.remove(paste0(data_dir, "seurat/sept_cornea_viral.h5Seurat"))  # Delete the old file

##########################################

## Viral Genes EDA
sept_cornea <- SeuratDisk::LoadH5Seurat(paste0(data_dir, "seurat/sept_cornea_viral_annot.h5Seurat"))

hsv_genes <- rownames(sept_cornea)[str_starts(rownames(sept_cornea), "HSV-")]

dge.kos_v_re <- FindMarkers(sept_cornea, 
                            features = hsv_genes,
                            ident.1 = "re", 
                            ident.2 = "kos", 
                            group.by = "condition", 
                            recorrect_umi = FALSE)

EnhancedVolcano::EnhancedVolcano(dge.kos_v_re, 
                                 lab = rownames(dge.kos_v_re), 
                                 #selectLab = hsv_genes,
                                 x = "avg_log2FC", y = "p_val_adj", 
                                 FCcutoff = 0.5, pCutoff =  1e-05, 
                                 title = "KOS vs RE",
                                 subtitle = "Sept Cornea All Cells",
                                 legendLabels = c("NS", expression(Log[2] ~ FC), "P-val", expression(P - value ~ and
                                                                                                     ~ Log[2] ~ FC)),
                                 legendPosition = "bottom",
                                 legendLabSize = 10,
                                 legendIconSize = 3,
                                 legendDropLevels = TRUE,
                                 colAlpha = 1, pointSize = 3, 
                                 drawConnectors = TRUE, widthConnectors = 0.5,
                                 labFace = "plain", boxedLabels = TRUE)

# Myeloid only
sept_cornea.mye = subset(sept_cornea, subset = "Mye" == cell_L2)

hsv_genes <- rownames(sept_cornea.mye)[str_starts(rownames(sept_cornea.mye), "HSV-")]

dge.kos_v_re.mye <- FindMarkers(sept_cornea.mye, 
                                features = hsv_genes,
                                ident.1 = "re", 
                                ident.2 = "kos", 
                                group.by = "condition", 
                                recorrect_umi = FALSE)

# NK only
sept_cornea.nk = subset(sept_cornea, subset = "NK" == cell_L2)

hsv_genes <- rownames(sept_cornea.nk)[str_starts(rownames(sept_cornea.nk), "HSV-")]

dge.kos_v_re.nk <- FindMarkers(sept_cornea.nk, 
                                features = hsv_genes,
                                ident.1 = "re", 
                                ident.2 = "kos", 
                                group.by = "condition", 
                                recorrect_umi = FALSE)

## Epithelial only
sept_cornea.epi = subset(sept_cornea, subset = "Epi" == cell_L2)
hsv_genes <- rownames(sept_cornea.epi)[str_starts(rownames(sept_cornea.epi), "HSV-")]

dge.kos_v_re.epi <- FindMarkers(sept_cornea.epi, 
                                features = hsv_genes,
                                ident.1 = "re", 
                                ident.2 = "kos", 
                                group.by = "condition", 
                                recorrect_umi = FALSE)

## Quick Violin Plot to look at the cell specificity of our viral DEGs
marker_genes_mye <- c("Itgam", "Adgre1", "Itgax",
                      "Itgam", # Cd11b - Myeloid marker (mouse)
                      "Adgre1", "Cd83", # F4/80, 
                      "Itgax", # DCs
                      "Cd14", "Cd68", "Fcgr1",  # Macs - note that Cd16 never comes up 
                      "Ly6c1", 
                      #"Cd74", # MHC-II mouse marker (used by Renthal 2022 to identify immune cells in TG)
                      "Ptgs2", "Irf5", "Nos2",  # Mouse M1 Mac Markers 
                      # "Stat1", "Retnla",  # Mouse M1 Mac Markers (less helpful)
                      #"Il12a", "Il23a", "Cd163",  # M1 vs M2 (M1: IL-12 and IL23 high with CD163 neg and M2 the opposite)
                      "Cd163",  # M2
                      #"Arg1", # M2a
                      "Socs3", "Cd86", # M2b
                      "Ccr2", "Slamf6",   #M2c
                      # "Tlr1", "Tlr8", "Scarb1", #M2c (less helpful)
                      "Vegfa",    # M2d, 
                      "Cx3cr1"  # Tissue-res Mac
)

VlnPlot(sept_cornea, assay = "SCT", 
        cols = c("black", "magenta"),
        features = hsv_genes, 
        group.by = "cell_L4", 
        split.by = "condition", split.plot = TRUE,
        stack = TRUE, flip = TRUE) +
  ggtitle("Sept Cornea")

#ggsave("viral_violin.svg", plot = last_plot(), path = figures_dir, width = 10, height = 14, units = "in", dpi = 300)
ggsave("viral_violin.png", plot = last_plot(), path = figures_dir, width = 7, height = 7, units = "in", dpi = 300, bg = "white")

VlnPlot(sept_cornea, assay = "SCT", 
        cols = c("black", "magenta"),
        features = c("HSV-RS1", "HSV-US1", "HSV-US10", "HSV-UL35"), 
        group.by = "cell_L4", 
        split.by = "condition", split.plot = TRUE,
        stack = TRUE, flip = TRUE) +
  ggtitle("Combined TG")

#ggsave("viral_violin.svg", plot = last_plot(), path = figures_dir, width = 10, height = 14, units = "in", dpi = 300)
ggsave("viral_violin2_date.png", plot = last_plot(), path = figures_dir, width = 7, height = 7, units = "in", dpi = 300, bg = "white")
