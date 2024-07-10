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

source("OmicsToolbox/marker_genes.R")

##########################################

# Read in Seurat Objects
seurat_obj <- SeuratDisk::LoadH5Seurat(paste0(data_dir, "sept_cornea.h5Seurat"))

## Factorizations
seurat_obj$condition <- factor(seurat_obj$condition, levels = c("control", "scratch", "kos", "re"))
seurat_obj$id <- factor(seurat_obj$id)
seurat_obj$date <- factor(seurat_obj$date)
seurat_obj$location <- factor(seurat_obj$location)
seurat_obj$celltype <- factor(seurat_obj$celltype)

##########################################


## Sept Cornea
VlnPlot(seurat_obj, assay = "SCT", 
        features = marker_genes_L2, 
        group.by = "seurat_clusters",
        #cols = scCustomize::DiscretePalette_scCustomize(length(marker_genes_L2), palette = "glasbey", shuffle = FALSE),
        stack = TRUE, flip = TRUE) + NoLegend() + 
  ggtitle("Sept Cornea (SCT)")

sept_cornea_L2_avgs <- AverageExpression(seurat_obj, assays = "SCT", layer = "scale.data", group.by = "seurat_clusters",
                                         features = c("Ptprc", "Epcam", "Itgam", "Cd74", "Cd68", "Cd3g", "Pax5", "Klrb1c", "Klrk1"))[[1]] %>% 
  as.matrix() %>% t() %>% 
  as_tibble(rownames = "cluster") %>% 
  mutate(designation = case_when(Epcam > 0.4 ~ "Epi", 
                                 Cd3g > 1 ~ "T", 
                                 (Itgam > 0) | ((Cd74 > 0) & (Pax5 < 1))  ~ "Mye", 
                                 Klrb1c > 0 ~ "NK", 
                                 Pax5 > 1 ~ "B", 
                                 (Ptprc > 0) | (Cd74> 0) ~ "other_imm",
                                 .default = "other"))

print(sept_cornea_L2_avgs, n = Inf)

seurat_obj$cell_L2 <- factor(sept_cornea_L2_avgs$designation[as.numeric(seurat_obj$seurat_clusters)], 
                              levels = c("Mye", "NK", "B", "T", "other_imm", "Epi", "other"))

seurat_obj$cell_L1 <- 
  factor(case_when(
    "Mye" == seurat_obj$cell_L2 ~ "immune", 
    "T" == seurat_obj$cell_L2 ~ "immune", 
    "B" == seurat_obj$cell_L2 ~ "immune", 
    "NK" == seurat_obj$cell_L2 ~ "immune", 
    "other_imm" == seurat_obj$cell_L2 ~ "immune", 
    "Epi" == seurat_obj$cell_L2 ~ "epithelial", 
    .default = "other"
  ), levels = c("immune", "epithelial", "other"))

## Cell L3 annotation
sept_cornea_cell_L4_mappings <- 
  seurat_obj@meta.data %>% 
  select(cell_L2, seurat_clusters) %>% 
  arrange(cell_L2, seurat_clusters) %>% 
  group_by(cell_L2) %>% 
  distinct() %>% 
  mutate(idx = row_number() - 1, 
         cell_L4 = paste0(cell_L2, "_", idx)) %>% 
  select(-idx) %>% 
  right_join(seurat_obj@meta.data %>% select(cell_L2, seurat_clusters) %>% rownames_to_column("barcode"), 
             join_by("cell_L2" == "cell_L2", "seurat_clusters" == "seurat_clusters")) %>% 
  column_to_rownames("barcode") %>% 
  select(cell_L4)

seurat_obj <- AddMetaData(seurat_obj, sept_cornea_cell_L4_mappings)
seurat_obj$cell_L4 <- factor(seurat_obj$cell_L4, 
                             levels =  str_sort(unique(seurat_obj$cell_L4), numeric = TRUE))


##########################################

## Write the new annotated version and delete our unannotated version
seurat_obj %>% SeuratDisk::SaveH5Seurat(paste0(data_dir, "sept_cornea_annot.h5Seurat"), overwrite = TRUE)

file.remove(paste0(data_dir, "sept_cornea.h5Seurat"))  # Delete the old file




