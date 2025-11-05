#Viral Counts Annotation and EDA

# Initialization ----
library(tidyverse)
library(ggpubr)
library(ggprism)
library(Seurat)
conflicted::conflict_prefer_all("dplyr", quiet = TRUE)

# DO NOT USE SEURAT V5 -> V5 objects do not play nicely yet with CellOracle in Python
options(Seurat.object.assay.version = "v3")

config_dir <- "./config/"
data_dir <- "../data_objects/"
figures_dir <- "../figures/"

## Parallell Processing ----
if(parallel::detectCores(logical=FALSE) > 3) {
  library(doParallel)
  
  num_cores <- parallel::detectCores(logical=FALSE)
  cl <- makePSOCKcluster(num_cores - 2)
  registerDoParallel(cl)
}

source("OmicsToolbox/marker_genes.R")
source("OmicsToolbox/seurat_norm_integrate.R")


# 10x objects -> Seurat Pipeline ----

configs <- c("combined_tg_viral.yaml")

for (config in configs) {
  loaded_config <- yaml::read_yaml(paste0(config_dir, config))

  seurat_norm_integrate(loaded_config, seurat_verbose = TRUE)
}



# Seurat Object Annotation ----

# Read in Seurat Object
combined_tg <- SeuratDisk::LoadH5Seurat(paste0(data_dir, "seurat/combined_tg_viral.h5Seurat"))

## Factorizations
combined_tg$condition <- factor(combined_tg$condition, levels = c("kos", "re"))
combined_tg$id <- factor(combined_tg$id)
combined_tg$date <- factor(combined_tg$date)
combined_tg$location <- factor(combined_tg$location)
combined_tg$celltype <- factor(combined_tg$celltype)


DimPlot(combined_tg, group.by = "seurat_clusters", 
        reduction = "umap", 
        label = TRUE, label.box = TRUE, label.size = 5, repel = FALSE, 
        cols = DiscretePalette(length(unique(combined_tg$seurat_clusters)), palette = "glasbey")
)

ggsave(paste0(figures_dir, "combined_tg_umap_seurat_clusters.png"), height = 7, width = 7)

VlnPlot(combined_tg, assay = "SCT",
        features = marker_genes_tg, 
        group.by = "seurat_clusters",
        stack = TRUE, flip = TRUE) + NoLegend() +
  ggtitle("Combined TG")

# B
FeaturePlot(combined_tg, 
            features = c("Pax5", "Cd19", "Cd20", "Ighg1", "Ighm", "Ighd"), 
            pt.size = 1.4)

# T
FeaturePlot(combined_tg, 
            features = c("Cd3d", "Cd3e"), 
            pt.size = 1.4)

## Crudely looking at this UMAP and marker genes, we have the following:
# B = 9
# T = 11
# NK = 5
# ? - 16 19 (remove these two)
# Myeloid for everything else

combined_tg <- subset(combined_tg, subset = seurat_clusters %in% c(1:15, 17, 18))

combined_tg_avgs <- AverageExpression(combined_tg, assays = "SCT", layer = "scale.data", group.by = "seurat_clusters",
                                      features = c("Cd3d", "Klrb1c", "Pax5", "Cd68"))[[1]] %>%
  as.matrix() %>% t() %>%
  as_tibble(rownames = "cluster") %>%
  mutate(designation = case_when(Pax5 > 0.3 ~ "B",
                                 Cd3d > 1 ~ "T",
                                 Klrb1c > 1 ~ "NK",
                                 #Cd68 > 0 ~ "Mye",
                                 .default = "Mye"
                                 #.default = cluster
  ),
  .after = cluster) %>% 
  mutate(cluster = sub("^g", "", cluster))

print(combined_tg_avgs, n = Inf)

combined_tg$cell_L1 <- factor("immune")

combined_tg$cell_L2 <- factor(combined_tg_avgs$designation[match(combined_tg$seurat_clusters, combined_tg_avgs$cluster)], 
                              levels = c("Mye", "NK", "B", "T"))


# Cell_L4 (e.g., "Mac_0", "Mac_1")
combined_tg_cell_L4_mappings <- 
  combined_tg@meta.data %>% 
  select(cell_L2, seurat_clusters) %>% 
  arrange(cell_L2, seurat_clusters) %>% 
  group_by(cell_L2) %>% 
  distinct() %>% 
  mutate(idx = row_number() - 1, 
         cell_L4 = paste0(cell_L2, "_", idx)) %>% 
  select(-idx) %>% 
  right_join(combined_tg@meta.data %>% select(cell_L2, seurat_clusters) %>% rownames_to_column("barcode"), 
             join_by("cell_L2" == "cell_L2", "seurat_clusters" == "seurat_clusters")) %>% 
  column_to_rownames("barcode") %>% 
  select(cell_L4)

combined_tg <- AddMetaData(combined_tg, combined_tg_cell_L4_mappings)
combined_tg$cell_L4 <- factor(combined_tg$cell_L4, 
                              levels =  str_sort(unique(combined_tg$cell_L4), numeric = TRUE))

## Take a quick view at the violin plot now to see if things line up as they should
VlnPlot(combined_tg, assay = "SCT",
        features = marker_genes_tg, 
        group.by = "cell_L4",
        stack = TRUE, flip = TRUE) + NoLegend() +
  ggtitle("Combined TG")


## Write the new annotated version and delete our unannotated version
combined_tg %>% SeuratDisk::SaveH5Seurat(paste0(data_dir, "combined_tg_viral_annot.h5Seurat"), overwrite = TRUE)

#file.remove(paste0(data_dir, "combined_tg_viral.h5Seurat"))  # Delete the old file


# Viral Genes EDA ----
combined_tg <- SeuratDisk::LoadH5Seurat(paste0(data_dir, "seurat/combined_tg_viral_annot.h5Seurat"))

hsv_genes <- rownames(combined_tg)[str_starts(rownames(combined_tg), "HSV-")]

dge.kos_v_re <- FindMarkers(combined_tg, 
                            features = hsv_genes,
                            ident.1 = "re", 
                            ident.2 = "kos", 
                            group.by = "condition", 
                            recorrect_umi = FALSE)

# Myeloid only
combined_tg.mye = subset(combined_tg, subset = "Mye" == cell_L2)
hsv_genes <- rownames(combined_tg.mye)[str_starts(rownames(combined_tg.mye), "HSV-")]

dge.kos_v_re.mye <- FindMarkers(combined_tg.mye, 
                                features = hsv_genes,
                                ident.1 = "re", 
                                ident.2 = "kos", 
                                group.by = "condition", 
                                recorrect_umi = FALSE)

EnhancedVolcano::EnhancedVolcano(dge.kos_v_re.mye, 
                                 lab = rownames(dge.kos_v_re.mye), 
                                 #selectLab = hsv_genes,
                                 x = "avg_log2FC", y = "p_val_adj", 
                                 FCcutoff = 0.5, pCutoff =  1e-05, 
                                 title = "KOS vs RE",
                                 subtitle = "TG Myeloid",
                                 legendLabels = c("NS", expression(Log[2] ~ FC), "P-val", expression(P - value ~ and
                                                                                                     ~ Log[2] ~ FC)),
                                 legendPosition = "bottom",
                                 legendLabSize = 10,
                                 legendIconSize = 3,
                                 legendDropLevels = TRUE,
                                 colAlpha = 1, pointSize = 3, 
                                 drawConnectors = TRUE, widthConnectors = 0.5,
                                 labFace = "plain", boxedLabels = TRUE)

#ggsave("viral_volcano.svg", plot = last_plot(), path = figures_dir, width = 12, height = 12, units = "in", dpi = 300)
ggsave("viral_volcano_mye.png", plot = last_plot(), path = figures_dir, width = 12, height = 12, units = "in", dpi = 300)

# NK only
combined_tg.nk = subset(combined_tg, subset = "NK" == cell_L2)
hsv_genes <- rownames(combined_tg.nk)[str_starts(rownames(combined_tg.nk), "HSV-")]

dge.kos_v_re.nk <- FindMarkers(combined_tg.nk, 
                                features = hsv_genes,
                                ident.1 = "re", 
                                ident.2 = "kos", 
                                group.by = "condition", 
                                recorrect_umi = FALSE)

EnhancedVolcano::EnhancedVolcano(dge.kos_v_re.nk, 
                                 lab = rownames(dge.kos_v_re.nk), 
                                 #selectLab = hsv_genes,
                                 x = "avg_log2FC", y = "p_val_adj", 
                                 FCcutoff = 0.5, pCutoff =  1e-05, 
                                 title = "KOS vs RE",
                                 subtitle = "TG NK Cells",
                                 legendLabels = c("NS", expression(Log[2] ~ FC), "P-val", expression(P - value ~ and
                                                                                                     ~ Log[2] ~ FC)),
                                 legendPosition = "bottom",
                                 legendLabSize = 10,
                                 legendIconSize = 3,
                                 legendDropLevels = TRUE,
                                 colAlpha = 1, pointSize = 3, 
                                 drawConnectors = TRUE, widthConnectors = 0.5,
                                 labFace = "plain", boxedLabels = TRUE)

#ggsave("viral_volcano.svg", plot = last_plot(), path = figures_dir, width = 12, height = 12, units = "in", dpi = 300)
ggsave("viral_volcano_nk.png", plot = last_plot(), path = figures_dir, width = 12, height = 12, units = "in", dpi = 300)

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

VlnPlot(combined_tg, assay = "SCT", 
        cols = c("black", "magenta"),
        features = c(rownames(dge.kos_v_re.mye), marker_genes_mye), 
        group.by = "cell_L4", 
        split.by = "condition", split.plot = TRUE,
        stack = TRUE, flip = TRUE) +
  ggtitle("Combined TG")

#ggsave("viral_violin.svg", plot = last_plot(), path = figures_dir, width = 10, height = 14, units = "in", dpi = 300)
ggsave("viral_violin.png", plot = last_plot(), path = figures_dir, width = 7, height = 7, units = "in", dpi = 300, bg = "white")

VlnPlot(combined_tg, assay = "SCT", 
        cols = c("black", "magenta"),
        features = c("HSV-RS1", "HSV-US1", "HSV-US10", "HSV-UL35", marker_genes_mye), 
        group.by = "cell_L4", 
        split.by = "condition", split.plot = TRUE,
        stack = TRUE, flip = TRUE) +
  ggtitle("Combined TG")

#ggsave("viral_violin.svg", plot = last_plot(), path = figures_dir, width = 10, height = 14, units = "in", dpi = 300)
ggsave("viral_violin2_date.png", plot = last_plot(), path = figures_dir, width = 7, height = 7, units = "in", dpi = 300, bg = "white")
