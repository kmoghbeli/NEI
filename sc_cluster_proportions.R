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

#######################################################################################################################################

dataset <- "sept_cornea"  # ** Most Cells - Use for Analysis
#dataset <- "combined_tg"

# Load the dataset
seurat_data <- SeuratDisk::LoadH5Seurat(paste0(data_dir, dataset, "_annot.h5Seurat"))  

# Sanity check for all of these (should all equal to 100)
# temp2 %>% group_by(condition) %>% summarise(summm = sum(prop))

# Cell_L2 by condition (just for all immune cells)
seurat_data@meta.data %>% 
  filter("immune" == cell_L1) %>% 
  group_by(date, condition, cell_L2) %>% 
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(prop = n / sum(n) * 100) %>%
  ggplot(aes(cell_L2, prop, fill = condition)) + 
  scale_fill_viridis_d() + 
  geom_col(position=position_dodge()) + 
  ggtitle(dataset) + 
  ggprism::theme_prism()

ggsave(paste0(figures_dir, dataset, "_props_immune_L2.pdf"), height = 7, width = 10)

# Cell_L4 by condition
seurat_data@meta.data %>% 
  filter("immune" == cell_L1, 
         condition %in% c("kos", "re")) %>% 
  group_by(date, condition, cell_L4) %>% 
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(prop = n / sum(n) * 100) %>%
  ggplot(aes(cell_L4, prop, fill = condition)) + 
  scale_fill_viridis_d(end = 0.75) + 
  geom_col(width = 0.55, position=position_dodge()) + 
  ggtitle(dataset) + 
  ggprism::theme_prism() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(paste0(figures_dir, dataset, "_props_immune_L4_kos_re.pdf"), height = 7, width = 10)


## Heatmap of DEGs between KOS and RE
# seurat_data.mye <- subset(seurat_data, subset = cell_L2 == "Mye" & condition %in% c("kos", "re"))
# kos_re_degs <- presto::wilcoxauc(seurat_data.mye, group_by = "condition", assay = "data", seurat_assay = "SCT") %>% 
#   filter("re" == group) %>% arrange(padj)
# 
# DoHeatmap(seurat_data.mye, features = kos_re_degs$feature[1:50], group.by = "condition")


## Violin Plot of Cell L3
VlnPlot(seurat_data, assay = "SCT", slot = "data", 
        features = marker_genes_L2, 
        group.by = "cell_L3",
        scCustomize::DiscretePalette_scCustomize(length(marker_genes_L2), palette = "varibow", shuffle = FALSE),
        stack = TRUE, flip = TRUE) + NoLegend() + 
  ggtitle("Sept Cornea (SCT)")

ggsave(paste0(figures_dir, dataset, "_markers_L3.pdf"), height = 10, width = 14)

## Umap Plot of Cell L3
DimPlot(seurat_data, 
        reduction = "umap.harmony", 
        cols = scCustomize::DiscretePalette_scCustomize(26, palette = "glasbey", shuffle = FALSE),
        group.by = c("cell_L3"))

ggsave(paste0(figures_dir, dataset, "_umap_L3.pdf"), height = 7, width = 10)

DimPlot(seurat_data, 
        reduction = "umap.harmony", 
        #cols = scCustomize::DiscretePalette_scCustomize(26, palette = "glasbey", shuffle = FALSE),
        group.by = c("condition"))

# Not a very helpful plot
#ggsave(paste0(figures_dir, dataset, "_umap_L3_split.pdf"), height = 7, width = 10)

## This code below gives you how many cells are in each cluster and which cell_L3 they map to
seurat_data@meta.data %>%
  select(seurat_clusters, condition) %>%
  table() %>%
  as.data.frame() %>% 
  pivot_wider(names_from = "condition", values_from = "Freq") %>% 
  left_join(seurat_data@meta.data %>% select(seurat_clusters, cell_L3) %>% distinct(seurat_clusters, cell_L3),
            join_by("seurat_clusters" == "seurat_clusters")) %>%
  select(seurat_clusters, cell_L3, everything()) %>% 
  print(n = Inf)
    


## Myeloid Cluster Plots
seurat_data.subset <- subset(seurat_data, subset = "Mye" == cell_L2 & condition %in% c("kos", "re") & nFeature_RNA > 500)

## Find the DEGs to keep

PrepSCTFindMarkers(seurat_data.subset)

seurat_data.subset$cell_L3 <- factor(seurat_data.subset$cell_L3, 
                                     levels = paste0("Mye_", seq(0,8)))
Idents(seurat_data.subset) <- seurat_data.subset$cell_L3
myeloid_degs <- FindAllMarkers(seurat_data.subset, only.pos = TRUE, assay = "SCT")
top_myeloid_degs <- myeloid_degs %>% 
  arrange(cluster) %>% 
  group_by(cluster) %>% 
  filter(p_val_adj <= 0.1) %>%     # Max FDR cut-off of 0.1
  arrange(p_val_adj) %>% 
  slice_head(n = 10)

mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256)
DoHeatmap(seurat_data.subset, 
          assay = "SCT", 
          raster = FALSE,
          features = top_myeloid_degs$gene, 
          size = 3) + 
  scale_fill_gradientn(colours = rev(mapal))

ggsave(paste0(figures_dir, dataset, "_top_myeloid_heatmap.pdf"), height = 10, width = 10)
# 
# immune_props_by_condition_cornea %>% 
#   ggplot(aes(cell_L3, prop, fill = condition)) + 
#   scale_fill_viridis_d() + 
#   geom_col(position=position_dodge(width = 0.75)) + 
#   theme_bw() + 
#   facet_grid(~ date) + 
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
# 
# ggsave(paste0(figures_dir, "props_immune_L3_cornea.pdf"), height = 7, width = 10)
# 
# 
# # Cell_L3 level - TG
# immune_props_by_condition_tg <- ni.sc@meta.data %>% 
#   filter("immune" == cell_L1, "tg" == location) %>% 
#   group_by(date, condition, cell_L3) %>% 
#   summarise(n = n(), .groups = "drop_last") %>%
#   mutate(prop = n / sum(n) * 100) 
# 
# immune_props_by_condition_tg %>% 
#   ggplot(aes(cell_L3, prop, fill = condition)) + 
#   scale_fill_viridis_d() + 
#   geom_col(position=position_dodge(width = 0.75)) + 
#   theme_bw() + 
#   facet_grid(~ date) + 
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
# 
# ggsave(paste0(figures_dir, "props_immune_L3_tg.pdf"), height = 7, width = 10)
# 
# 
# # Cell_L3 level - TG vs Cornea
# immune_props_by_condition_location <- ni.sc@meta.data %>% 
#   filter("immune" == cell_L1) %>% 
#   group_by(location, condition, cell_L3) %>% 
#   summarise(n = n(), .groups = "drop_last") %>%
#   mutate(prop = n / sum(n) * 100) 
# 
# immune_props_by_condition_location %>% 
#   ggplot(aes(cell_L3, prop, fill = condition)) + 
#   scale_fill_viridis_d() + 
#   geom_col(position=position_dodge(width = 0.75)) + 
#   theme_bw() + 
#   facet_grid(~ location) + 
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
# 
# ggsave(paste0(figures_dir, "props_immune_L3_location.pdf"), height = 7, width = 10)
# 
# 
# # immune_props_by_date_and_condition <- ni.sc@meta.data %>% 
# #   filter("immune" == cell_L1, "neuron_date" != date) %>% 
# #   group_by(date, condition, cell_L3) %>% 
# #   summarise(n = n(), .groups = "drop_last") %>%
# #   mutate(prop = n / sum(n) * 100) 
# # 
# # immune_props_by_date_and_condition %>% 
# #   ggplot(aes(cell_L3, prop, fill = condition)) + 
# #   scale_fill_viridis_d() + 
# #   geom_col(position=position_dodge(width = 0.75)) + 
# #   facet_grid(~ date) + 
# #   theme_bw()
# # 
# # immune_props_by_date_and_condition %>% 
# #   filter(date %in% c("20230413", "20230706"), condition %in% c("control", "scratch", "re")) %>% 
# #   ggplot(aes(cell_L3, prop, fill = condition)) + 
# #   scale_fill_viridis_d() + 
# #   geom_col(position=position_dodge(width = 0.75)) + 
# #   facet_grid(~ date) + 
# #   theme_bw()
# 
# 
# # Similarly for the epithelial cell clusters:
#   
# epi_props_by_condition <- ni.sc@meta.data %>% 
#   filter("epithelial" == cell_L1, "cornea" == location) %>% 
#   group_by(date, condition, cell_L3) %>% 
#   summarise(n = n(), .groups = "drop_last") %>%
#   mutate(prop = n / sum(n) * 100) 
# 
# epi_props_by_condition %>% 
#   ggplot(aes(cell_L3, prop, fill = condition)) + 
#   scale_fill_viridis_d() + 
#   geom_col(position=position_dodge()) + 
#   facet_grid(~ date) + 
#   theme_bw()
# 
# ggsave(paste0(figures_dir, "props_epi_L3.pdf"))#, height = 7, width = 10)