library(tidyverse)
library(ggpubr)
library(ggprism)
library(Seurat)
conflicted::conflict_prefer_all("dplyr", quiet = TRUE)

# DO NOT USE SEURAT V5 -> V5 objects do not play nicely yet with CellOracle in Python
options(Seurat.object.assay.version = "v3")


kaveh_colors1 <- c("#a6611a", "#dfc27d", "#f5f5f5", "#80cdc1", "#018571")
kaveh_colors2 <- c("#40004b", "#762a83", "#9970ab", "#c2a5cf", "#e7d4e8", "#f7f7f7", "#d9f0d3", "#a6dba0", "#5aae61", "#1b7837", "#00441b")
kaveh_colors3 <- c("#ff0000", "#f5f5f5","#0000ff")


config_dir <- "./config/"
data_dir <- "./data/"
model_dir <- "./models/"
results_dir <- "./results/"
figures_dir <- "./figures/"

# Parallell Processing
if(parallel::detectCores(logical=FALSE) > 3) {
  library(doParallel)
  
  num_cores <- parallel::detectCores(logical=FALSE)
  cl <- makePSOCKcluster(num_cores - 2)
  registerDoParallel(cl)
}

marker_genes_L1 <- c("Epcam", "Ptprc")

# Source: https://www.cellsignal.com/pathways/immune-cell-markers-mouse
# https://www.thermofisher.com/us/en/home/life-science/cell-analysis/cell-analysis-learning-center/immunology-at-work/macrophage-cell-overview.html
marker_genes_L2 <- c("Epcam", "Ptprc", 
                     "Itgam", # Cd11b - Myeloid marker (mouse)
                     "Adgre1", # F4/80, 
                     "Cd74", # MHC-II mouse marker (used by Renthal 2022 to identify immune cells in TG)
                     "Pax5",  # B 
                     "Ighd", "Cd27", # Naive (mouse) B cell markers (IgD+, CD27-)
                     "Cd14", "Cd68",  # Macs - note that Cd16 never comes up 
                     "Cox2", "Irf5", "Nos2", "Stat1",  # Mouse M1 Mac Markers
                     "Cd3d", "Cd3e", "Cd3g", "Cd4", "Cd8a",   # T
                     "Cd44", "Sell", # Naive (mouse) T cell markers (CD44 -, CD62L +)
                     "Fcgr3", "Il17a", "Ifng",    # delta-gamma T
                     "Klrb1c", "Klrk1", "Gzma", "Gzmb", "Prf1", # NK
                     "Itga2", "Ncam1",  #NK-T
                     "Pdgfra",  # Fibroblasts
                     "Vwf")     # Endothelial

# CD64 - mouse cornea macrophages (maybe activation marker, i.e., low in M0)
# CCR2 maybe also an activation marker
# Ly6c - inflammatory monocytes
# Cd68
# Nature Communications: Qie et al ("integrated proteomic... macrophages") paper and supplemental data: macrophage.mousprotein.cn


##########################################

# Read in Seurat Objects
## July Cornea
july_cornea <- SeuratDisk::LoadH5Seurat(paste0(data_dir, "july_cornea_sct.h5Seurat"))
sept_cornea <- SeuratDisk::LoadH5Seurat(paste0(data_dir, "sept_cornea_sct.h5Seurat"))  # ** Most Cells - Use for Analysis
june_tg <- SeuratDisk::LoadH5Seurat(paste0(data_dir, "june_tg_sct.h5Seurat"))
aug_tg <- SeuratDisk::LoadH5Seurat(paste0(data_dir, "aug_tg_sct.h5Seurat"))


# july_cornea.marker_data_wide <- GetAssayData(july_cornea, layer = "scale.data", assay="SCT") %>% 
#   as.matrix() %>% t() %>% as.data.frame() %>% rownames_to_column("barcode") %>% 
#   select(barcode, any_of(marker_genes_L2)) 
# 
# july_cornea.marker_data_long <- july_cornea.marker_data_wide %>% 
#   pivot_longer(cols = !any_of(c("barcode")), names_to = "gene", values_to = "expr")


##########################################

# Clustree Analysis 
# Note that this is not really needed after it is done once for a dataset. 
#  For our individual datasets above, a cluster resolution of 1.0 seems to be fine

library(clustree)  # note that this call is needed otherwise the calls to "clustree" function will error (see: https://github.com/lazappi/clustree/issues/14)
  
clustree(july_cornea, prefix = "SCT_snn_res.", 
                   node_colour = "sc3_stability", 
                   show_axis = TRUE) + ggtitle("July Cornea (SCT)")

clustree(sept_cornea, prefix = "SCT_snn_res.", 
                   node_colour = "sc3_stability", 
                   show_axis = TRUE) + ggtitle("Sept Cornea (SCT)")

clustree(june_tg, prefix = "SCT_snn_res.", 
                   node_colour = "sc3_stability", 
                   show_axis = TRUE) + ggtitle("June TG (SCT)")

clustree(aug_tg, prefix = "SCT_snn_res.", 
                   node_colour = "sc3_stability", 
                   show_axis = TRUE) + ggtitle("Aug TG (SCT)")

# Cluster resolution of ~1.0 seems to be ok
july_cornea <- FindClusters(july_cornea, resolution = 1.0)
sept_cornea <- FindClusters(sept_cornea, resolution = 1.0)
june_tg <- FindClusters(june_tg, resolution = 1.0)
aug_tg <- FindClusters(aug_tg, resolution = 1.0)

##########################################
# L2 Markers

## July Cornea
VlnPlot(july_cornea, assay = "SCT", 
        features = marker_genes_L2, 
        group.by = "seurat_clusters",
        scCustomize::DiscretePalette_scCustomize(length(marker_genes_L2), palette = "glasbey", shuffle = FALSE),
        stack = TRUE, flip = TRUE) + NoLegend() + 
  ggtitle("July Cornea (SCT)")

july_cornea_L2_avgs <- AverageExpression(july_cornea, assays = "SCT", layer = "scale.data", group.by = "seurat_clusters",
                                         features = c("Ptprc", "Epcam", "Itgam", "Cd74", "Cd3g", "Pax5", "Klrb1c", "Klrk1"))[[1]] %>% 
  as.matrix() %>% t() %>% 
  as_tibble(rownames = "cluster") %>% 
  mutate(designation = case_when(Epcam > 0 ~ "Epi", 
                                 Cd3g > 0 ~ "T", 
                                 (Itgam > 0) & (Cd74 > 0) ~ "Mye", 
                                 Klrb1c > 0 ~ "NK", 
                                 #Pax5 > 0 ~ "B", 
                                 (Ptprc > 0) | (Cd74> 0) ~ "Immune",
                                 .default = "other"))

print(july_cornea_L2_avgs, n = Inf)

july_cornea$cell_L2 <- factor(july_cornea_L2_avgs$designation[as.numeric(july_cornea$seurat_clusters)], 
                                    levels = c("Mye", "T", "B", "NK", "Immune", "Epi", "other"))

july_cornea$cell_L1 <- 
  factor(case_when(
    "Mye" == july_cornea$cell_L2 ~ "immune", 
    "T" == july_cornea$cell_L2 ~ "immune", 
    "B" == july_cornea$cell_L2 ~ "immune", 
    "NK" == july_cornea$cell_L2 ~ "immune", 
    "Immune" == july_cornea$cell_L2 ~ "immune", 
    "Epi" == july_cornea$cell_L2 ~ "epithelial", 
    .default = "other"
  ), levels = c("immune", "epithelial", "other"))

## Sept Cornea
VlnPlot(sept_cornea, assay = "SCT", 
        features = marker_genes_L2, 
        group.by = "seurat_clusters",
        scCustomize::DiscretePalette_scCustomize(length(marker_genes_L2), palette = "glasbey", shuffle = FALSE),
        stack = TRUE, flip = TRUE) + NoLegend() + 
  ggtitle("Sept Cornea (SCT)")

sept_cornea_L2_avgs <- AverageExpression(sept_cornea, assays = "SCT", layer = "scale.data", group.by = "seurat_clusters",
                                         features = c("Ptprc", "Epcam", "Itgam", "Cd74", "Cd68", "Cd3g", "Pax5", "Klrb1c", "Klrk1"))[[1]] %>% 
  as.matrix() %>% t() %>% 
  as_tibble(rownames = "cluster") %>% 
  mutate(designation = case_when(Epcam > 0 ~ "Epi", 
                                 Cd3g > 0 ~ "T", 
                                 Itgam > 0 ~ "Mye", 
                                 Klrb1c > 0 ~ "NK", 
                                 Pax5 > 0 ~ "B", 
                                 (Ptprc > 0) | (Cd74> 0) ~ "Immune",
                                 .default = "other"))

print(sept_cornea_L2_avgs, n = Inf)

sept_cornea$cell_L2 <- factor(sept_cornea_L2_avgs$designation[as.numeric(sept_cornea$seurat_clusters)], 
                              levels = c("Mye", "T", "B", "NK", "Immune", "Epi", "other"))

sept_cornea$cell_L1 <- 
  factor(case_when(
    "Mye" == sept_cornea$cell_L2 ~ "immune", 
    "T" == sept_cornea$cell_L2 ~ "immune", 
    "B" == sept_cornea$cell_L2 ~ "immune", 
    "NK" == sept_cornea$cell_L2 ~ "immune", 
    "Immune" == sept_cornea$cell_L2 ~ "immune", 
    "Epi" == sept_cornea$cell_L2 ~ "epithelial", 
    .default = "other"
  ), levels = c("immune", "epithelial", "other"))

## Cell L3 annotation
sept_cornea_cell_L3_mappings <- 
  sept_cornea@meta.data %>% 
  select(cell_L2, seurat_clusters) %>% 
  arrange(cell_L2, seurat_clusters) %>% 
  group_by(cell_L2) %>% 
  distinct() %>% 
  mutate(idx = row_number() - 1, 
         cell_L3 = paste0(cell_L2, "_", idx)) %>% 
  select(-idx) %>% 
  right_join(sept_cornea@meta.data %>% select(cell_L2, seurat_clusters) %>% rownames_to_column("barcode"), 
             join_by("cell_L2" == "cell_L2", "seurat_clusters" == "seurat_clusters")) %>% 
  column_to_rownames("barcode") %>% 
  select(cell_L3)

sept_cornea <- AddMetaData(sept_cornea, sept_cornea_cell_L3_mappings)



## June TG
VlnPlot(june_tg, assay = "SCT", 
        features = marker_genes_L2, 
        group.by = "seurat_clusters",
        scCustomize::DiscretePalette_scCustomize(length(marker_genes_L2), palette = "glasbey", shuffle = FALSE),
        stack = TRUE, flip = TRUE) + NoLegend() + 
  ggtitle("June TG (SCT)")

june_tg_L2_avgs <- AverageExpression(june_tg, assays = "SCT", layer = "scale.data", group.by = "seurat_clusters",
                                    features = c("Ptprc", "Cd74", "Cd3g", "Pax5", "Klrb1c", "Klrk1"))[[1]] %>% 
  as.matrix() %>% t() %>% 
  as_tibble(rownames = "cluster") %>% 
  mutate(designation = case_when(Cd3g > 0 ~ "T", 
                                 Klrb1c > 0 ~ "NK", 
                                 Pax5 > 0 ~ "B", 
                                 .default = "Immune"))

print(june_tg_L2_avgs, n = Inf)

june_tg$cell_L2 <- factor(june_tg_L2_avgs$designation[as.numeric(june_tg$seurat_clusters)], 
                         levels = c("Mye", "T", "B", "NK", "Immune", "Epi", "other"))

june_tg$cell_L1 <- factor("immune")


## Aug TG
VlnPlot(aug_tg, assay = "SCT", 
        features = marker_genes_L2, 
        group.by = "seurat_clusters",
        scCustomize::DiscretePalette_scCustomize(length(marker_genes_L2), palette = "glasbey", shuffle = FALSE),
        stack = TRUE, flip = TRUE) + NoLegend() + 
  ggtitle("Aug TG (SCT)")

aug_tg_L2_avgs <- AverageExpression(aug_tg, assays = "SCT", layer = "scale.data", group.by = "seurat_clusters",
                                         features = c("Cd74", "Cd3g", "Pax5", "Klrb1c", "Klrk1"))[[1]] %>% 
  as.matrix() %>% t() %>% 
  as_tibble(rownames = "cluster") %>% 
  mutate(designation = case_when(Pax5 > 0 ~ "B", 
                                 Cd3g > 2 ~ "T", 
                                 Klrk1 > 0 ~ "NK", 
                                 .default = "Immune"))

print(aug_tg_L2_avgs, n = Inf)

aug_tg$cell_L2 <- factor(aug_tg_L2_avgs$designation[as.numeric(aug_tg$seurat_clusters)], 
                              levels = c("Mye", "T", "B", "NK", "Immune", "Epi", "other"))

aug_tg$cell_L1 <-  factor("immune")

##########################################

# Factorizations and write to disk

## Factorizations
july_cornea$condition <- factor(july_cornea$condition, levels = c("control", "scratch", "kos", "re"))
july_cornea$id <- factor(july_cornea$id)
july_cornea$date <- factor(july_cornea$date)
july_cornea$location <- factor(july_cornea$location)
july_cornea$celltype <- factor(july_cornea$celltype)

sept_cornea$condition <- factor(sept_cornea$condition, levels = c("control", "scratch", "kos", "re"))
sept_cornea$id <- factor(sept_cornea$id)
sept_cornea$date <- factor(sept_cornea$date)
sept_cornea$location <- factor(sept_cornea$location)
sept_cornea$celltype <- factor(sept_cornea$celltype)

june_tg$condition <- factor(june_tg$condition, levels = c("control", "scratch", "kos", "re"))
june_tg$id <- factor(june_tg$id)
june_tg$date <- factor(june_tg$date)
june_tg$location <- factor(june_tg$location)
june_tg$celltype <- factor(june_tg$celltype)

aug_tg$condition <- factor(aug_tg$condition, levels = c("control", "scratch", "kos", "re"))
aug_tg$id <- factor(aug_tg$id)
aug_tg$date <- factor(aug_tg$date)
aug_tg$location <- factor(aug_tg$location)
aug_tg$celltype <- factor(aug_tg$celltype)


## Write the new annotated version and delete our unannotated version
july_cornea %>% SeuratDisk::SaveH5Seurat(paste0(data_dir, "july_cornea_sct.h5Seurat"), overwrite = TRUE)
sept_cornea %>% SeuratDisk::SaveH5Seurat(paste0(data_dir, "sept_cornea_sct.h5Seurat"), overwrite = TRUE)
june_tg %>% SeuratDisk::SaveH5Seurat(paste0(data_dir, "june_tg_sct.h5Seurat"), overwrite = TRUE)
aug_tg %>% SeuratDisk::SaveH5Seurat(paste0(data_dir, "aug_tg_sct.h5Seurat"), overwrite = TRUE)



