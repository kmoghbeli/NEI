library(tidyverse)
library(ggpubr)
library(ggprism)
library(Seurat)
library(harmony)
library(leiden)

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

##########################################

neuron_counts <- readr::read_csv("../../../_DATASETS/NeuroImmune/neuron_counts_all_kallisto_tximport_lengthscaledtpm.csv",
                                 show_col_types = FALSE) %>%
  column_to_rownames(var = "gene")

colnames(neuron_counts) <- sub("neuron_", "", colnames(neuron_counts))

neuron_counts <- neuron_counts[, !is.na(colSums(neuron_counts)) & colSums(neuron_counts) > 0] %>%
  as.matrix() %>%
  Matrix::Matrix(sparse = T)

neuron_metadata <- readr::read_csv("../../../_DATASETS/NeuroImmune/neuron_metadata_with_filenames.csv",
                                   show_col_types = FALSE) %>%
  select(-filename) %>% 
  column_to_rownames(var = "cell")

neuron_obj <- CreateSeuratObject(neuron_counts,
                                 project = "NeuroImmune", 
                                 meta.data = neuron_metadata)

## Update/Set some Metadata
neuron_obj$condition <- factor(neuron_obj$condition, levels = c("control", "scratch", "kos", "re"))
neuron_obj$location <- "corneal_afferents"
neuron_obj$celltype <- "neuron"
neuron_obj$celltypeloc <- "neuron"
neuron_obj$timepoint <- ""
neuron_obj$cell_L1 <- neuron_obj$cell_L2 <- neuron_obj$cell_L3 <- neuron_obj$cell_L4 <- "neuron"

neuron_obj <- SCTransform(neuron_obj, 
                          vst.flavor = "v2",
                          #vars.to.regress = "percent.mt",
                          return.only.var.genes = FALSE)

neuron_obj <- RunPCA(neuron_obj, assay = "SCT", npcs = 50)

neuron_obj <- RunHarmony(neuron_obj,
                         assay.use = "SCT",
                         reduction = "pca",
                         dims.use = 1:50, 
                         lambda = NULL,  # Ridge regression penalty - when set to NULL, harmony tries to estimate
                         kmeans_init_nstart=20, 
                         kmeans_init_iter_max=100, 
                         group.by.vars = c("condition", "batch_folder"), 
                         reduction.save = "harmony", 
                         plot_convergence = TRUE)

neuron_obj <- RunUMAP(neuron_obj, assay = "SCT", reduction = "harmony", dims = 1:50)
neuron_obj <- FindNeighbors(neuron_obj, assay = "SCT", reduction = "harmony", dims = 1:50)

neuron_obj <- FindClusters(neuron_obj, resolution = 1.0, method = "igraph", algorithm = "Leiden")

neuron_obj <- PrepSCTFindMarkers(neuron_obj)

# ## quick plots to show that batch correction worked
# p1 <- DimPlot(neuron_obj, group.by = "seurat_clusters") + theme(legend.position = "bottom")
# p2 <- DimPlot(neuron_obj, group.by = "condition") + theme(legend.position = "bottom")
# p3 <- DimPlot(neuron_obj, group.by = "batch_folder") + theme(legend.position = "bottom") + guides(colour = guide_legend(ncol = 2))
# p1 + p2 + p3

neuron_obj %>% SeuratDisk::SaveH5Seurat(paste0(data_dir, "neurons_all_conditions.h5Seurat"), overwrite = TRUE)

