# Seurat Pipeline

library(tidyverse)
library(ggpubr)
library(ggprism)
library(Seurat)
conflicted::conflict_prefer_all("dplyr", quiet = TRUE)

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

# DO NOT USE SEURAT V5 -> V5 objects do not play nicely yet with CellOracle in Python
options(Seurat.object.assay.version = "v3")

## MARKER DEFINITIONS

# https://www.cellsignal.com/pathways/immune-cell-markers-mouse
# Renthal 2022: 
#   - "Sparc" (non-neuronal marker gene
#   - Cd74 = immune
#   - non-neuronal subtype marker genes (Satglia = Apoe, Schwann cells = Mpz, fibroblasts = Dcn or Mgp, immune cells = Cd74, and vascular = Igfbp7)

marker_genes_L1 <- c("Epcam", "Ptprc")

marker_genes_L2 <- c("Epcam", "Ptprc", 
                     "Itgam", # Cd11b - Myeloid marker (mouse)
                     "Pax5",  # B 
                     "Cd14", "Cd16", "Cd68", # Mac
                     "Cd3d", "Cd3e", "Cd3g", "Cd4", "Cd8a", # T
                     "Klrb1c", "Klrk1", "Gzma", "Gzmb", "Prf1", # NK
                     "Pdgfra",  # Fibroblasts
                     "Vwf")     # Endothelial

config <- yaml::read_yaml(paste0(config_dir, "sept_cornea_sct.yaml"))


##### CLUSTER IDENTITY EXPLORATION
##### THIS PART HAS TO BE CUSTOMIZED FOR EACH CONFIG/RUN ####### 

# Read in Seurat Object
harmonized_seurat <- SeuratDisk::LoadH5Seurat(paste0(data_dir, config$filename, ".h5Seurat"))

### Level 1 ("cell_L1") designations (Immune, Epithelial, etc)
##### THRESHOLDS OF 0.7 (Epcam) and 1.5 (Ptprc) are for July Cornea only

# Visualize
VlnPlot(harmonized_seurat, assay = "RNA", 
        features = marker_genes_L1, 
        group.by = "seurat_clusters",
        scCustomize::DiscretePalette_scCustomize(length(marker_genes_L1), palette = "glasbey", shuffle = FALSE),
        stack = TRUE, flip = TRUE) + NoLegend()

VlnPlot(harmonized_seurat, assay = "RNA", 
        features = marker_genes_L2, 
        group.by = "seurat_clusters",
        scCustomize::DiscretePalette_scCustomize(length(marker_genes_L2), palette = "glasbey", shuffle = FALSE),
        stack = TRUE, flip = TRUE) + NoLegend()


epcam_ptprc_avgs <- AverageExpression(harmonized_seurat, assays = "RNA", layer = "data", group.by = "seurat_clusters",
                                      features = c("Epcam", "Ptprc"))[[1]] %>% as.matrix() %>% t() %>% as_tibble(rownames = "cluster")

if ("sct" == config$normalization) {
  ## SCT thresholds
  epcam_ptprc_avgs <- epcam_ptprc_avgs%>% 
    mutate(designation = case_when(Epcam > 1 ~ "epithelial", 
                                   Ptprc > 2 ~ "immune", 
                                   .default = "other"))
} else {
  ## Log normalization thresholds
  
  epcam_ptprc_avgs <- epcam_ptprc_avgs%>% 
    mutate(designation = case_when(Epcam > 1 ~ "epithelial", 
                                   Ptprc > 2 ~ "immune", 
                                   .default = "other"))
}


harmonized_seurat$cell_L1 <- factor(epcam_ptprc_avgs$designation[as.numeric(harmonized_seurat$seurat_clusters)], 
                                    levels = c("immune", "epithelial", "other"))



### Level 2 ("cell_L2") designations (Myeloid, B, T, NK, etc)
  
  
VlnPlot(harmonized_seurat, assay = "RNA", 
        features = marker_genes_L2, 
        group.by = "seurat_clusters",
        scCustomize::DiscretePalette_scCustomize(length(marker_genes_L2), palette = "glasbey", shuffle = FALSE),
        stack = TRUE, flip = TRUE) + NoLegend()

ggsave(paste0(figures_dir, "violin_cell_L2.pdf"), height = 5, width = 5)

immune_L2_avgs <- AverageExpression(harmonized_seurat, assays = "RNA", layer = "data", group.by = "seurat_clusters",
                                      features = c("Epcam", "Itgam", "Cd3g", "Pax5", "Klrb1c", "Klrk1"))[[1]] %>% as.matrix() %>% t() %>% as_tibble(rownames = "cluster") %>% 
  mutate(designation = case_when(Epcam > 1 ~ "Epi", 
                                 Cd3g > 5 ~ "T", 
                                 Itgam > 1.5 ~ "Mye", 
                                 Klrb1c > 1 ~ "NK", 
                                 Klrk1 > 2 ~ "NKT", 
                                 Pax5 > 5 ~ "B", 
                                 .default = "other"))

harmonized_seurat$cell_L2 <- factor(immune_L2_avgs$designation[as.numeric(harmonized_seurat$seurat_clusters)], 
                                    levels = c("Mye", "T", "B", "NK", "NKT", "Epi", "other"))

# Sanity Check
harmonized_seurat@meta.data %>% select(cell_L1, cell_L2) %>% 
  group_by(cell_L1, cell_L2) %>% 
  distinct() %>% 
  arrange(cell_L1, cell_L2)

## Factorizations
harmonized_seurat$condition <- factor(harmonized_seurat$condition, 
                                      levels = c("control", "scratch", "kos", "re"))
harmonized_seurat$id <- factor(harmonized_seurat$id)
harmonized_seurat$date <- factor(harmonized_seurat$date)
harmonized_seurat$location <- factor(harmonized_seurat$location)
harmonized_seurat$celltype <- factor(harmonized_seurat$celltype)


## Write the new annotated version and delete our unannotated version
harmonized_seurat %>% SeuratDisk::SaveH5Seurat(paste0(data_dir, config$filename, ".h5Seurat"), overwrite = TRUE)

##########################################
## Cell Proportions

# Cell_L2 level
harmonized_seurat@meta.data %>% 
  filter("immune" == cell_L1) %>% 
  group_by(location, condition, cell_L2) %>% 
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(prop = n / sum(n) * 100) %>%
  ggplot(aes(cell_L2, prop, fill = condition)) + 
  scale_fill_viridis_d() + 
  geom_col(position=position_dodge()) + 
  #facet_grid(~ location) + 
  ggprism::theme_prism() + theme(axis.title.x = element_blank()) + 
  ylab("Percentage of All Immune Cells for Condition")

ggsave(paste0(figures_dir, "props_immune_L2_location.pdf"), height = 7, width = 7)

##########################################
# Cluster Analysis - double-negative T cells??

harmonized_seurat$clus_L2 <- paste0(harmonized_seurat$seurat_clusters, "_", harmonized_seurat$cell_L2)
harmonized_seurat$clus_L2 <- factor(harmonized_seurat$clus_L2, 
                                    levels = gtools::mixedsort(unique(harmonized_seurat$clus_L2)))

marker_data_wide <- GetAssayData(harmonized_seurat, layer = "data", assay="RNA") %>% 
  as.matrix() %>% t() %>% as.data.frame() %>% rownames_to_column("barcode") %>% 
  select(barcode, any_of(marker_genes_L2)) %>% 
  left_join(tibble(barcode = names(harmonized_seurat$clus_L2), cluster = harmonized_seurat$clus_L2), 
            by = c("barcode" = "barcode"))

marker_data_long <- marker_data_wide %>% 
  pivot_longer(cols = !any_of(c("barcode", "cluster")), names_to = "gene", values_to = "expr") %>% 
  filter(expr > 0)

marker_data_long %>% 
  filter(gene %in% c("Ptprc", "Epcam")) %>% 
  ggplot(aes(expr, color = gene)) + 
  #geom_histogram(aes(y=..density..), bins = 100, alpha = 0.25) + 
  geom_density(aes(y = ..count..)) + 
  #geom_density(aes(y = ..scaled..)) + 
  facet_wrap(~cluster, scales = "free_y") + 
  theme(legend.title = element_blank())

marker_data_long %>% 
  filter(gene %in% c("Cd3d", "Cd4", "Cd8a")) %>% 
  ggplot(aes(expr, color = gene)) + 
  #geom_histogram(aes(y=..density..), bins = 100, alpha = 0.25) + 
  geom_density(aes(y = ..count..)) + 
  #geom_density(aes(y = ..scaled..)) + 
  facet_wrap(~cluster, scales = "free_y") + 
  theme(legend.title = element_blank())


marker_data_wide %>% ggplot(aes(Ptprc, Epcam, color = cluster)) + geom_point()

marker_data_wide %>% 
  ggplot(aes(Cd3d, Cd4)) + 
  #geom_point(shape = ".") + 
  geom_point(alpha = 0.25) + 
  facet_wrap(~cluster)

marker_data_wide %>% ggplot(aes(Cd3d, Cd8a, color = cluster)) + geom_point()

VlnPlot(harmonized_seurat, assay = "RNA", 
        features = marker_genes_L2, 
        group.by = "seurat_clusters",
        scCustomize::DiscretePalette_scCustomize(length(marker_genes_L2), palette = "glasbey", shuffle = FALSE),
        stack = TRUE, flip = TRUE) + NoLegend()

##########################################
# 
# 
# marker_genes <- c("Epcam", "Ptprc", 
#                   #"Trpv1", "Ly6g", "Siglech", 
#                   "Bst2", "Cd83", 
#                   #"Itgax", 
#                   "Itgam", # Cd11b - Myeloid marker (mouse)
#                   "Ly6c1", 
#                   "Cd80", "Cd86", "Mrc1", "Nos2", "H2-Ab1",
#                   "Adgre1", "Cd163", "Cd14", # Macs
#                   "Pax5", "Cd19",   # B 
#                   #"Cd27",
#                   #"Cd19", "Cd20", "Ighd", "Cd27", "Msa41", "Pax5", # B
#                   "Cd3d", "Klrb1c", "Klrk1", "Gzma", "Gzmb", "Prf1", "Cd69",   # T/NK
#                   #"Tnfrsf17", "Sdc1",
#                   "Pdgfra", "Lum", 
#                   #"Kera", 
#                   "Vwf")
# 
# Idents(harmonized_seurat) <- harmonized_seurat$seurat_clusters
# VlnPlot(harmonized_seurat, assay = "RNA", 
#         features = marker_genes, 
#         scCustomize::DiscretePalette_scCustomize(length(marker_genes), palette = "glasbey", shuffle = FALSE),
#         stack = TRUE, flip = TRUE) + NoLegend()
# 
# ggsave(paste0(figures_dir, "violin_markers.pdf"), height = 7, width = 7)
# 
# # Idents(harmonized_seurat) <- harmonized_seurat$seurat_clusters
# # DotPlot(harmonized_seurat, assay = "RNA", 
# #         features = marker_genes) + coord_flip()
# # 
# # ggsave(paste0(figures_dir, "dotplot_markers.pdf"), height = 7, width = 10)
# 
# 
# harmonized_seurat$condition <- factor(harmonized_seurat$condition, 
#                                       levels = c("control", "scratch", "kos", "re"))
# 
# harmonized_seurat$cell_L1 <- 
#   case_match(harmonized_seurat$seurat_clusters, 
#              as.factor(c(0:12, 15:17, 19:21, 24:27)) ~ "immune", 
#              as.factor(c(13, 14, 18)) ~ "epithelial", 
#              as.factor(22) ~ "endothelial", 
#              as.factor(23) ~ "stromal/fibroblast")
# 
# harmonized_seurat$cell_L2 <- 
#   case_match(harmonized_seurat$seurat_clusters, 
#              as.factor(c(0, 1, 2, 3, 6, 8, 9, 10, 12, 15, 16, 17, 19, 20, 21, 24, 26)) ~ "Mac", 
#              as.factor(c(7, 11, 25)) ~ "B", 
#              as.factor(c(5)) ~ "T", 
#              as.factor(c(4, 27)) ~ "NK",
#              as.factor(c(13, 14, 18)) ~ "Epi", 
#              as.factor(c(22)) ~ "endothelial", 
#              as.factor(c(23)) ~ "stromal/fibroblast")
# 
# 
# harmonized_seurat$cell_L3 <- harmonized_seurat$seurat_clusters
# levels(harmonized_seurat$cell_L3) <- c("Mac1", "Mac2", "Mac3", "Mac4", "NK1", "T1", 
#                                        "Mac5", "B1", "Mac6", "Mac7", "Mac8", "B2", 
#                                        "Mac9", "Epi1", "Epi2", "Mac10", "Mac11", "Mac12", 
#                                        "Epi3", "Mac13", "Mac14", "Mac15", "Endo", "Strom", 
#                                        "Mac16", "B3", "Mac17", "NK2")
# 
# # Re-order the levels
# harmonized_seurat$cell_L2 <- factor(harmonized_seurat$cell_L2, 
#                                     levels = c("Mac", "T", "B", "NK", "Epi", "endothelial", "stromal/fibroblast"))
# 
# harmonized_seurat$cell_L3 <- factor(harmonized_seurat$cell_L3, 
#                                     levels = c(paste0("Mac", 1:17), 
#                                                "T1",   
#                                                "B1", "B2", "B3", 
#                                                "NK1", "NK2", 
#                                                "Epi1", "Epi2", "Epi3", 
#                                                "Endo", "Strom"))
# 
# # Sanity Check
# harmonized_seurat@meta.data %>% select(cell_L1, cell_L2, cell_L3) %>% 
#   group_by(cell_L1, cell_L2, cell_L3) %>% 
#   distinct() %>% 
#   arrange(cell_L1, cell_L2, cell_L3)
# 
# # Idents(harmonized_seurat) <- harmonized_seurat$cell_ident
# # VlnPlot(harmonized_seurat, assay = "RNA", 
# #         features = marker_genes, 
# #         scCustomize::DiscretePalette_scCustomize(length(marker_genes), palette = "varibow", shuffle = FALSE),
# #         stack = TRUE, flip = TRUE) + NoLegend()
# 
# # Keep only immune, epithelial, and neuron for downstream analyses
# Idents(harmonized_seurat) <- harmonized_seurat$cell_L1
# 
# ni.sc <- subset(harmonized_seurat, idents = c("immune", "epithelial"))

# 
# # Write to disk (and also write AnnData version)
# harmonized_seurat %>% readr::write_rds(paste0(data_dir, "neuroimmune_harmonized_seurat.rds"))
# ni.sc %>% readr::write_rds(paste0(data_dir, "ni_sc.rds"))
