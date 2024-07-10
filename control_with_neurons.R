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

# Read in Sept Cornea, Combined TG, and then Neurons
# sept_cornea <- SeuratDisk::LoadH5Seurat(paste0(data_dir, "sept_cornea_annot.h5Seurat"))
# sept_cornea_control <- subset(sept_cornea, subset = "control" == condition)
# sept_cornea_control$celltypeloc <- paste0("cornea_", sept_cornea_control$cell_L2)
# rm(sept_cornea) # free memory since these objects are so large

combined_tg <- SeuratDisk::LoadH5Seurat(paste0(data_dir, "combined_tg_annot.h5Seurat"))
combined_tg_control <- subset(combined_tg, subset = "control" == condition)
#combined_tg_control$celltypeloc <- paste0("tg_", combined_tg_control$cell_L2)
combined_tg_control$celltypeloc <- combined_tg_control$cell_L2
rm(combined_tg) # free memory since these objects are so large

##########################################
# Neuron Data - Read in SmartSeq count matrix (should have cells as columns and genes as rows)
neuron_counts <- readr::read_csv("../../../_DATASETS/NeuroImmune/control_corneal_afferent_counts_filtered.csv",
                                 show_col_types = FALSE) %>%
  column_to_rownames(var = "gene")

neuron_counts <- neuron_counts[, colSums(neuron_counts) > 0] %>%
  as.matrix() %>%
  Matrix::Matrix(sparse = T)

neuron_obj <- CreateSeuratObject(neuron_counts,
                                 project = "NeuroImmune")

message("Neurons: ", nrow(neuron_obj), " genes x ", ncol(neuron_obj), " cells")

neuron_obj$id <- "HandpickedControlCornealAfferents"
neuron_obj$date <- "neuron_date"
neuron_obj$condition <- "control"
neuron_obj$location <- "corneal_afferents"
neuron_obj$celltype <- "neuron"
neuron_obj$celltypeloc <- "neuron"
neuron_obj$timepoint <- ""

neuron_obj$cell_L1 <- neuron_obj$cell_L2 <- neuron_obj$cell_L3 <- neuron_obj$cell_L4 <- "neuron"

# Skipping QC on neurons since they were hand-picked single cells

neuron_obj <- SCTransform(neuron_obj, 
                          vst.flavor = "v2",
                          #vars.to.regress = "percent.mt",
                          return.only.var.genes = FALSE)

##########################################
# Merge everything together 
seurat_obj <- merge(x = neuron_obj,
                    y = c(combined_tg_control),
                    #y = c(sept_cornea_control, combined_tg_control),
                    merge.data = TRUE)

# seurat_obj <- RunPCA(object = seurat_obj, assay = "SCT", npcs = 50, verbose = seurat_verbose)
# seurat_obj <- RunUMAP(object = seurat_obj, assay = "SCT", reduction = "harmony", dims = 1:50)

# Do we need to run this?
seurat_obj <- PrepSCTFindMarkers(seurat_obj)

## Factorizations
seurat_obj$condition <- factor(seurat_obj$condition)
seurat_obj$id <- factor(seurat_obj$id)
seurat_obj$date <- factor(seurat_obj$date)
seurat_obj$location <- factor(seurat_obj$location)
seurat_obj$celltype <- factor(seurat_obj$celltype)
seurat_obj$cell_L1 <- factor(seurat_obj$cell_L1)
seurat_obj$cell_L2 <- factor(seurat_obj$cell_L2)
seurat_obj$cell_L3 <- factor(seurat_obj$cell_L3)
seurat_obj$cell_L4 <- factor(seurat_obj$cell_L4)
seurat_obj$celltypeloc <- factor(seurat_obj$celltypeloc)

##########################################d
VlnPlot(seurat_obj, assay = "SCT", 
        features = c(marker_genes_L2, "Trpv1"), 
        group.by = "celltypeloc",
        #cols = scCustomize::DiscretePalette_scCustomize(length(marker_genes_L2), palette = "glasbey", shuffle = FALSE),
        stack = TRUE, flip = TRUE) + NoLegend() + 
  ggtitle("Control with Neurons")

##########################################

## Write the new annotated version and delete our unannotated version
seurat_obj %>% SeuratDisk::SaveH5Seurat(paste0(data_dir, "tg_control_with_neurons.h5Seurat"), overwrite = TRUE)




