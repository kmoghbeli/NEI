# Seurat Pipeline

library(tidyverse)

source("OmicsToolbox/seurat_norm_integrate.R")

conflicted::conflict_prefer_all("dplyr", quiet = TRUE)

#reticulate::py_install("leidenalg")
#leidenalg <- reticulate::import("leidenalg")

kaveh_colors1 <- c("#a6611a", "#dfc27d", "#f5f5f5", "#80cdc1", "#018571")
kaveh_colors2 <- c("#40004b", "#762a83", "#9970ab", "#c2a5cf", "#e7d4e8", "#f7f7f7", "#d9f0d3", "#a6dba0", "#5aae61", "#1b7837", "#00441b")
kaveh_colors3 <- c("#ff0000", "#f5f5f5","#0000ff")


config_dir <- "./config/"
data_dir <- "./data/"
model_dir <- "./models/"
results_dir <- "./results/"
figures_dir <- "./figures/"


#configs <- c("june_tg_sct.yaml", "aug_tg_sct.yaml", "july_cornea_sct.yaml", "sept_cornea_sct.yaml")
#configs <- c("sept_cornea.yaml")
#configs <- c("control_sept_cornea_combined_tg.yaml")
configs <- c("combined_tg.yaml")
#configs <- c("combined_cornea.yaml")

for (config in configs) {
  loaded_config <- yaml::read_yaml(paste0(config_dir, config))
  
  seurat_norm_integrate(loaded_config, seurat_verbose = TRUE)
}
