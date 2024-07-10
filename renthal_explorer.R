library(tidyverse)
library(ggpubr)
library(ggprism)
library(Seurat)
library(leiden)
library(harmony)
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

DATASETS_dir <- "../../../_DATASETS/"

# Parallell Processing
if(parallel::detectCores(logical=FALSE) > 3) {
  library(doParallel)
  
  num_cores <- parallel::detectCores(logical=FALSE)
  cl <- makePSOCKcluster(num_cores - 2)
  registerDoParallel(cl)
}

source("OmicsToolbox/marker_genes.R")

##########################################

# Read in Seurat Object
renthal_seurat <- readr::read_rds(paste0(DATASETS_dir, "Renthal/Renthal_Seurat_v0.2.Rds"))

renthal_mouse <- subset(renthal_seurat, subset = "Mouse" == Species)

renthal_mouse@meta.data %>% 
  distinct(publication, platform, Species, Strain, Age, Sex) %>% 
  rownames_to_column("barcode") %>% 
  select(-barcode)

table(renthal_mouse$final_annotation) %>% as.data.frame() %>% rename("Annotation" = Var1)

DimPlot(renthal_mouse, 
        cols = DiscretePalette(length(unique(renthal_seurat$final_annotation)), palette = "glasbey", shuffle = FALSE), 
        group.by = "final_annotation")


##########################################
## Trigeminal Ganglion Mouse dataset
counts <- readr::read_rds(paste0(DATASETS_dir, "Renthal/tg_neuron_2022/snRNA-seq_mouse_raw_counts.RDS")) 

metadata <- readr::read_delim(paste0(DATASETS_dir, "Renthal/tg_neuron_2022/GSE197289_snRNA-seq_mouse_barcode_meta.csv"))

## Some of the cells have duplicates in the metadata
duplicates <- metadata$V1[which(duplicated(metadata$V1))]
metadata %>%
  filter(V1 %in% dups, class == "non-neuron") %>%
  select(V1, class, cellID, subtype, library) %>%
  group_by(V1) %>% 
  summarise(count = n()) %>% 
  arrange(V1) 

## Subset just the non-neuronal cells and remove the few with duplicate entries for the 2 libraries
metadata_non_neuronal <- metadata %>% 
  filter("non-neuron" == class, !(V1 %in% duplicates)) %>% 
  column_to_rownames("V1")
  
counts_non_neuronal <- counts[, rownames(metadata_non_neuronal)]

obj <- CreateSeuratObject(counts = counts_non_neuronal, meta.data = metadata_non_neuronal, min.cells = 3, min.features = 200)

# No need to do pre-filtering since it was done on this dataset already

seurat_objs <- SplitObject(obj, split.by = "sample_name")

for (i in 1:length(seurat_objs)) {
  seurat_objs[[i]] <- SCTransform(seurat_objs[[i]], 
                                  vst.flavor = "v2",
                                  vars.to.regress = "percent.mt",
                                  return.only.var.genes = FALSE,
                                  verbose = TRUE)
}

integ.features <- SelectIntegrationFeatures(object.list = seurat_objs, nfeatures = 3000)

## Merge
merged_seurat <- merge(x = seurat_objs[[1]],
                       y = seurat_objs[2:length(seurat_objs)],
                       merge.data = TRUE)

VariableFeatures(merged_seurat) <- integ.features

merged_seurat <- RunPCA(object = merged_seurat, assay = "SCT", npcs = 50, verbose = TRUE)

merged_seurat <- RunHarmony(object = merged_seurat,
                            assay.use = "SCT",
                            reduction = "pca",
                            dims.use = 1:50, 
                            lambda = NULL,  # Ridge regression penalty - when set to NULL, harmony tries to estimate
                            kmeans_init_nstart=20, kmeans_init_iter_max=100, 
                            group.by.vars = c("model", "sample_name"), 
                            reduction.save = "harmony", 
                            plot_convergence = TRUE)

merged_seurat <- RunUMAP(object = merged_seurat, assay = "SCT", reduction = "harmony", dims = 1:50)
merged_seurat <- FindNeighbors(object = merged_seurat, assay = "SCT", reduction = "harmony", dims = 1:50, verbose = TRUE)

merged_seurat <- FindClusters(object = merged_seurat, resolution = 1.0, method = "igraph", algorithm = "Leiden")

merged_seurat <- PrepSCTFindMarkers(object = merged_seurat, verbose = TRUE)

merged_seurat %>% readr::write_rds(paste0(data_dir, "renthal_tg_non_neuronal.rds"))

VlnPlot(merged_seurat, assay = "SCT", 
        features = marker_genes_L2, 
        group.by = "seurat_clusters",
        #cols = scCustomize::DiscretePalette_scCustomize(length(marker_genes_L2), palette = "glasbey", shuffle = FALSE),
        stack = TRUE, flip = TRUE) + NoLegend() + 
  ggtitle("Renthal Non-Neurons")

merged_seurat.naive_immune <- subset(merged_seurat, subset = seurat_clusters == 23 & model == "Naive")

merged_seurat.naive_immune <- RunUMAP(object = merged_seurat.naive_immune, assay = "SCT", reduction = "harmony", dims = 1:50)
merged_seurat.naive_immune <- FindNeighbors(object = merged_seurat.naive_immune, assay = "SCT", reduction = "harmony", dims = 1:50, verbose = TRUE)

merged_seurat.naive_immune <- FindClusters(object = merged_seurat.naive_immune, resolution = 2.0, method = "igraph", algorithm = "Leiden")

merged_seurat.naive_immune <- PrepSCTFindMarkers(object = merged_seurat.naive_immune, verbose = TRUE)

VlnPlot(merged_seurat.naive_immune, assay = "SCT", 
        features = c("Ptprc", "Itgam", "Adgre1", "Cd74", "Cd68", "Pax5", "Cd3d", "Klrb1c", "Klrk1"), 
        group.by = "seurat_clusters",
        stack = TRUE, flip = TRUE) + NoLegend() + 
  ggtitle("Renthal Non-Neurons")

##########################################


# Idents(combined_tg) <- combined_tg$cell_L2
# not_other <- levels(combined_tg$cell_L2)["other" != levels(combined_tg$cell_L2)]
# combined_tg_no_other <- subset(combined_tg, idents = not_other)
# 
# DimPlot(combined_tg_no_other, group.by = "cell_L2", reduction = "umap", label = TRUE, label.box = TRUE, label.size = 3, repel = TRUE) + 
#   ggtitle("Trigeminal Ganglion") + NoLegend()
# 
# ggsave(paste0(figures_dir, "combined_tg_umap.png"), height = 7, width = 7)


# Sanity check for all of these (should all equal to 100)
# combined_tg@meta.data %>% 
#   group_by(condition, cell_L2) %>% 
#   summarise(n = n(), .groups = "drop_last") %>%
#   mutate(prop = n / sum(n) * 100) %>% 
#   summarise(summm = sum(prop))

# Cell L2 and L3 by condition (just for all immune cells)
combined_tg@meta.data %>% 
  group_by(condition, cell_L2) %>% 
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(prop = n / sum(n) * 100) %>% 
  ggplot(aes(cell_L2, prop, fill = condition)) + 
  scale_fill_viridis_d() + 
  geom_col(position=position_dodge()) + 
  ggprism::theme_prism() + 
  labs(x = "", y = "Proportion of all immune cells", title = "Trigeminal Ganglion")

ggsave(paste0(figures_dir, "combined_tg_props_immune_L2.png"), height = 7, width = 14)

combined_tg@meta.data %>% 
  # filter("immune" == cell_L1, 
  #        condition %in% c("kos", "re")) %>% 
  group_by(date, condition, cell_L2) %>% 
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(prop = n / sum(n) * 100) %>%
  ggplot(aes(cell_L2, prop, fill = condition)) + 
  scale_fill_viridis_d() + 
  geom_col(width = 0.55, position=position_dodge()) + 
  facet_grid(~ date) + 
  ggprism::theme_prism() + 
  #theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  labs(x = "", y = "Proportion of all immune cells", title = "Trigeminal Ganglion")

ggsave(paste0(figures_dir, "combined_tg_props_immune_L2_date.png"), height = 7, width = 14)


combined_tg@meta.data %>% 
  group_by(condition, cell_L3) %>% 
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(prop = n / sum(n) * 100) %>% 
  ggplot(aes(cell_L3, prop, fill = condition)) + 
  scale_fill_viridis_d() + 
  geom_col(position=position_dodge()) + 
  ggprism::theme_prism() + 
  labs(x = "", y = "Proportion of all immune cells", title = "Trigeminal Ganglion")

ggsave(paste0(figures_dir, "combined_tg_props_immune_L3.png"), height = 7, width = 14)

combined_tg@meta.data %>% 
  # filter("immune" == cell_L1, 
  #        condition %in% c("kos", "re")) %>% 
  group_by(date, condition, cell_L3) %>% 
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(prop = n / sum(n) * 100) %>%
  ggplot(aes(cell_L3, prop, fill = condition)) + 
  scale_fill_viridis_d() + 
  geom_col(width = 0.55, position=position_dodge()) + 
  facet_grid(~ date) + 
  ggprism::theme_prism() + 
  #theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  labs(x = "", y = "Proportion of all immune cells", title = "Trigeminal Ganglion")

ggsave(paste0(figures_dir, "combined_tg_props_immune_L3_date.png"), height = 7, width = 20)


combined_tg@meta.data %>% 
  filter("immune" == cell_L1, 
          condition %in% c("kos", "re")) %>% 
  group_by(date, condition, cell_L3) %>% 
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(prop = n / sum(n) * 100) %>%
  ggplot(aes(cell_L3, prop, fill = condition)) + 
  scale_fill_viridis_d(end = 0.75) + 
  geom_col(width = 0.55, position=position_dodge()) + 
  facet_grid(~ date) + 
  ggprism::theme_prism() + 
  #theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  labs(x = "", y = "Proportion of all immune cells", title = "Trigeminal Ganglion")

ggsave(paste0(figures_dir, "combined_tg_props_immune_L3_date_re_kos.png"), height = 7, width = 21)


table(combined_tg$seurat_clusters, combined_tg$condition)

##########################################

## Write the new annotated version and delete our unannotated version
combined_tg %>% SeuratDisk::SaveH5Seurat(paste0(data_dir, "combined_tg.h5Seurat"), overwrite = TRUE)



