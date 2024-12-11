## NeuroImmune SLIDE

library(tidyverse)
library(Seurat)
#conflicted::conflict_prefer_all("dplyr", quiet = TRUE)

kaveh_colors1 <- c("#a6611a", "#dfc27d", "#f5f5f5", "#80cdc1", "#018571")
kaveh_colors2 <- c("#40004b", "#762a83", "#9970ab", "#c2a5cf", "#e7d4e8", "#f7f7f7", "#d9f0d3", "#a6dba0", "#5aae61", "#1b7837", "#00441b")
kaveh_colors3 <- c("#ff0000", "#f5f5f5","#0000ff")


config_dir <- "./config/"
data_dir <- "./data/"
model_dir <- "./models/"
results_dir <- "./results/"
figures_dir <- "./figures/"

# Parallell Processing
num_cores <- future::availableCores()
if(num_cores > 3) {
  library(doParallel)
  
  cl <- makePSOCKcluster(num_cores - 2)
  registerDoParallel(cl)
}

source("slide_helper_functions.R")

# Run Slide Immune analyses on Sept Cornea, June TG, Aug TG
# - all immune, myeloid, and T cells, epithelial

# Datasets to prep
# slide_comparisons <- expand_grid(dataset = c("july_cornea_sct", "sept_cornea_sct", "june_tg_sct", "aug_tg_sct"), 
#                                  metadata_key = c("cell_L1", "cell_L2"), 
#                                  metadata_val = c("immune", "epithelial", "Mye", "T", "B"), 
#                                  compare1 = c("control", "scratch", "kos", "re"), 
#                                  compare2 = c("control", "scratch", "kos", "re")) %>% 
#   filter(compare1 != compare2, 
#          (grepl("cornea", dataset) & (metadata_val %in% c("immune", "epithelial", "Mye", "T"))) | 
#            (grepl("tg", dataset) & (metadata_val %in% c("immune", "B", "T"))), 
#          ("cell_L1" == metadata_key & metadata_val %in% c("immune", "epithelial")) | 
#            ("cell_L2" == metadata_key & metadata_val %in% c("Mye", "T", "B"))) %>% 
#   rowwise() %>% 
#   mutate(id = paste(dataset, metadata_key, metadata_val, min(compare1, compare2), max(compare1, compare2), sep = "_")) %>% 
#   distinct(id, .keep_all = TRUE) %>% 
#   select(-id)

slide_comparisons <- expand_grid(dataset = c("combined_tg"), 
                                 metadata_key = c("cell_L2"), 
                                 metadata_val = list(c("Mac", "T"), c("T")), 
                                 compare1 = c("kos"), 
                                 compare2 = c("re")) %>% 
  rowwise() %>% 
  mutate(id = paste(dataset, metadata_key, paste(metadata_val, sep = "_", collapse = "_"), min(compare1, compare2), max(compare1, compare2), sep = "_")) %>% 
  distinct(id, .keep_all = TRUE) %>% 
  select(-id)


### PIPELINE ### 
for (curr_dataset in unique(slide_comparisons$dataset)) {
  
  dataset_comps <- slide_comparisons %>% filter(dataset == curr_dataset)

  config <- yaml::read_yaml(paste0(config_dir, curr_dataset, ".yaml"))
  
  seurat_obj <- SeuratDisk::LoadH5Seurat(paste0(data_dir, config$filename, ".h5Seurat"))
    
  for (i in 1:nrow(dataset_comps)) { 
    comp <- dataset_comps[i, ]
    
    id <- paste0(curr_dataset, "_", paste(unlist(comp$metadata_val), collapse = "_"), "_", toupper(comp$compare1), "v", toupper(comp$compare2))
    
    print(paste0("Prepping: ", id))
    
    data_subset <- subsetNIData(seurat_obj, comp$metadata_key, comp$metadata_val[[1]], c(comp$compare1, comp$compare2))
    
    prepSLIDE(data = data_subset, 
              out_path = paste0("./slide_runs/", id, "/"),
              lambda = c(0.1, 0.5, 1.0),   # higher = sparser LFs (i.e., more loadings are zero) (0.5 - 1.0)
              delta =  c(0.01, 0.1),  # higher = fewer LFs returned from the unsupervised part (0.01, 0.1)
              spec = 0.3,     # higher = fewer significant LFs returned
              thresh_fdr = 0.2)
  }
}


#SLIDE::checkDataParams(yaml::yaml.load_file("./slide_runs/combined_tg_Mac_T_KOSvRE/er.yaml"))

