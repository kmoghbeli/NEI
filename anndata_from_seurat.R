library(tidyverse)
library(ggpubr)
library(ggprism)
library(Seurat)
conflicted::conflict_prefer_all("dplyr", quiet = TRUE)

config_dir <- "./config/"
data_dir <- "../data_objects/"

# Parallell Processing
if(parallel::detectCores(logical=FALSE) > 3) {
  library(doParallel)
  
  num_cores <- parallel::detectCores(logical=FALSE)
  cl <- makePSOCKcluster(num_cores - 2)
  registerDoParallel(cl)
}

#######################################################################################################################################

dataset <- "combined_tg_viral_annot"

## 1) Load the full annotated, normalized dataset

seurat_data <- SeuratDisk::LoadH5Seurat(paste0(data_dir, "seurat/", dataset, ".h5Seurat"))  

## 2) Make the Cell L2 annotations (e.g., Mye, B, NK, etc) into strings (not factors)
seurat_data$cell_L2 <- as.character(seurat_data$cell_L2)

## 3) Do some fudging of count storage because of how "Convert" will store stuff in the AnnData object
##  - store the "data" slot of "SCT" into "scale.data"
##  - store the raw counts (from "RNA" assay "counts") into "SCT" "data" so that Convert will store it into "raw.X"
seurat_data[["SCT"]]$scale.data <- seurat_data[["SCT"]]$data %>% as.matrix()
seurat_data[["SCT"]]$data <- seurat_data[["RNA"]]$counts

seurat_data %>% SeuratDisk::SaveH5Seurat(paste0(data_dir, dataset, "_temp_for_celloracle.h5Seurat"), overwrite = TRUE)
SeuratDisk::Convert(source = paste0(data_dir, dataset, "_temp_for_celloracle.h5Seurat"), 
                    dest = paste0(data_dir, "anndata/", dataset, ".h5ad"), 
                    overwrite = TRUE)
file.remove(paste0(data_dir, dataset, "_temp_for_celloracle.h5Seurat"))  # Delete the temp file
