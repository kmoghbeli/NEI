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

#dataset <- "sept_cornea_sct"  # ** Most Cells - Use for Analysis for Cornea
dataset <- "combined_tg"
conditions <- c("kos", "re")

## 1) Load the full annotated, normalized dataset

seurat_data <- SeuratDisk::LoadH5Seurat(paste0(data_dir, dataset, ".h5Seurat"))  


## 2) Subset just our cells of interest amongst the KOS vs RE conditions

#seurat_data.subset <- subset(seurat_data, subset = "Mye" == cell_L2 & condition %in% conditions & nFeature_RNA > 500)
seurat_data.subset <- subset(seurat_data, subset = "Mac" == cell_L2 & condition %in% conditions & nFeature_RNA > 500)


## 3) find the DEGs to keep

seurat_data.subset <- PrepSCTFindMarkers(seurat_data.subset)

degs <- FindMarkers(seurat_data.subset, group.by = "condition", ident.1 = "re", ident.2 = "kos") %>% 
  filter(p_val_adj <= 0.1) %>%     # Max FDR cut-off of 0.1
  arrange(p_val_adj) %>% 
  slice_head(n = 3000) %>%     # Keep no more than top 3k genes
  rownames_to_column("gene")


## 4) Load the SLIDE results for this dataset subset
base_slide_path <- "slide_runs/combined_tg_Mac_KOSvRE_D0.01_L0.5_S0.3/"
lfs <- c(87, 27, 139, 18, 147, 31, 107, 157)

lf_features <- readr::read_rds("slide_runs/combined_tg_Mac_KOSvRE_D0.01_L0.5_S0.3/plotSigGenes_data.RDS") %>%
  mutate(names = ifelse(grepl("^X.+Rik$", names), sub("^X", "", names), names)) %>%
  mutate(names = gsub("\\.", "-", names)) %>%  # For some reason ER/SLIDE results swap out "-" for "." in our gene names 
  rename(gene = names) %>% 
  filter(lf_num %in% lfs)
  
## 5) Downsample the number of cells (needed to run CellOracle if dataset is too large e.g., greater than 3k cells)
re_cells = Cells(seurat_data.subset)[which(seurat_data.subset$condition == "re")]
kos_cells = Cells(seurat_data.subset)[which(seurat_data.subset$condition == "kos")]

downsampled_re_cells = sample(re_cells, size = 1500)
downsampled_kos_cells = sample(kos_cells, size = 1500)


## X) Final Subset
seurat_data.final_subset <- 
  subset(seurat_data,
         subset = "Mac" == cell_L2 & condition %in% c("kos", "re") & nFeature_RNA > 500, 
         cells = c(downsampled_re_cells, downsampled_kos_cells), 
         features = c(degs$gene, lf_features$gene))


## 2) Make the Cell L3 annotations into strings (not factors)
seurat_data.final_subset$cell_L3 <- as.character(seurat_data.final_subset$cell_L3)

## 2a) Add some annotations for SLIDE LF genes and Hocomo TFs
hocomo_tfs <- readr::read_tsv("hocomocoV12_TFs.txt") %>% pull(var = 1)

feature_metadata <- tibble(gene = rownames(seurat_data.final_subset[["SCT"]])) %>% 
  mutate(is_tf = str_to_lower(gene) %in% str_to_lower(hocomo_tfs), 
         is_slide_lf = gene %in% lf_features$gene) %>% 
  column_to_rownames("gene")

seurat_data.final_subset[["SCT"]] <- AddMetaData(seurat_data.final_subset[["SCT"]], metadata = feature_metadata)


## 3) Do some fudging of count storage because of how "Convert" will store stuff in the AnnData object
##  - store the "data" slot of "SCT" into "scale.data"
##  - store the raw counts (from "RNA" assay "counts") into "SCT" "data" so that Convert will store it into "raw.X"
seurat_data.final_subset[["SCT"]]$scale.data <- seurat_data.final_subset[["SCT"]]$data %>% as.matrix()
seurat_data.final_subset[["SCT"]]$counts <- seurat_data.final_subset[["RNA"]]$counts

seurat_data.final_subset %>% SeuratDisk::SaveH5Seurat(paste0(data_dir, dataset, "_temp_for_celloracle.h5Seurat"), overwrite = TRUE)
SeuratDisk::Convert(source = paste0(data_dir, dataset, "_temp_for_celloracle.h5Seurat"), 
                    dest = paste0(data_dir, "../../NEICellOracle/", dataset, "_mac_kos_vs_re.h5ad"), 
                    overwrite = TRUE)
file.remove(paste0(data_dir, dataset, "_temp_for_celloracle.h5Seurat"))  # Delete the temp file




## Code to read HoCoMoCo V12 TFs
# lines <- readLines("H12CORE-MOUSE_annotation.jsonl")
# lines <- lapply(lines, jsonlite::fromJSON)
# lines <- lapply(lines, unlist)
# tfs <- bind_rows(lines)
# tfs %>% select(original_motif.tf) %>% distinct(original_motif.tf) %>% readr::write_tsv("hocomocoV12_TFs.txt", col_names = FALSE)



### Some Plots to Explore the Data
# DimPlot(seurat_data.subset, 
#         reduction = "umap", 
#         cols = scCustomize::DiscretePalette_scCustomize(9, palette = "glasbey", shuffle = FALSE),
#         group.by = c("cell_L3"), 
#         split.by = c("condition"))
# 
# ggsave(paste0(figures_dir, dataset, "_umap_L3_split.pdf"), height = 7, width = 10)
# 
# seurat_data@meta.data %>% 
#   filter(condition %in% conditions) %>%
#   group_by(condition, cell_L3) %>% 
#   summarise(n = n(), .groups = "drop_last") %>%
#   mutate(prop = n / sum(n) * 100) %>% 
#   ggplot(aes(cell_L3, prop, fill = condition)) + 
#   scale_fill_viridis_d(begin = 0.14, end = 0.75) + 
#   geom_col(position=position_dodge()) + 
#   ggpubr::theme_pubr() + 
#   labs(x = "", y = "Proportion of all immune cells", title = "Trigeminal Ganglion")
# 
# ggsave(paste0(figures_dir, "combined_tg_props_immune_L3_kos_re.png"), height = 7, width = 14)
# 
# seurat_data@meta.data %>% 
#   filter(condition %in% conditions) %>%
#   group_by(condition, seurat_clusters) %>% 
#   summarise(n = n(), .groups = "drop_last") %>%
#   mutate(prop = n / sum(n) * 100) %>% 
#   ggplot(aes(seurat_clusters, prop, fill = condition)) + 
#   scale_fill_viridis_d(begin = 0.14, end = 0.75) + 
#   geom_col(position=position_dodge()) + 
#   ggpubr::theme_pubr() + 
#   labs(x = "", y = "Proportion of all immune cells", title = "Trigeminal Ganglion")
# 
# ggsave(paste0(figures_dir, "combined_tg_props_immune_seurat_clusters_kos_re.png"), height = 7, width = 14)