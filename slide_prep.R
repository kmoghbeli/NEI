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

### HELPER FUNCTIONS ### 
#' Reduce sparsity of the data matrix by requiring samples to have a minimun number of features 
#' and features to be expressed in a minimum number of samples
#'
#' @param data the sparse data matrix (assuming rows are features and columns are samples)
#' @param sample_thresh samples with fewer than sample_thresh number of features will be filtered out
#' @param feature_thresh features expressed in fewer than feature_thresh number of samples will be filtered out
#' @export
sparsityFiltering <- function(data, sample_thresh, feature_thresh){
  print(paste0("Original dataframe dimension is ", nrow(data), " by ", ncol(data)))
  
  binary_matrix <- data != 0
  
  #filter out features expressed in fewer than feature_thresh number of samples
  i <- rowSums(binary_matrix, na.rm=TRUE) >= feature_thresh
  
  #filtered out samples with fewer than sample_thresh number of features
  j <- colSums(binary_matrix, na.rm=TRUE) >= sample_thresh
  
  filtered <- data[i ,j]
  
  print(paste0("Filtered dataframe dimension is ", nrow(filtered), " by ", ncol(filtered)))
  
  return(filtered)
}

#subsetNIData(seurat_obj, comp$metadata_key, comp$metadata_val[[1]], c(comp$compare1, comp$compare2))

subsetNIData <- function(seurat_obj, 
                         celltype_metadata_key, 
                         celltype_metadata_val, 
                         conditions, 
                         exclude_mito_ribo_genes = TRUE) {
  
  max_cells_per_condition <- 10000000
  max_features <- 5000
  
  print(paste0("subsetNIData: Initial dataset - ", ncol(seurat_obj), " cells and ", nrow(seurat_obj), " features"))
  
  ## 1) Remove mitochondrial and ribosomal genes
  if (TRUE == exclude_mito_ribo_genes) {
    genes_to_keep <- which(!grepl("^MT-|^RPS|^RPL|^mRP|^ATP", 
                                  Features(seurat_obj), 
                                  ignore.case = TRUE))
    
    seurat_obj <- subset(seurat_obj, features = Features(seurat_obj)[genes_to_keep])
  }  
  
  print(paste0("subsetNIData: Removed mito/ribo genes - ", ncol(seurat_obj), " cells and ", nrow(seurat_obj), " features"))
  
  ## 2) First Zero Filter
  data <- GetAssayData(seurat_obj, layer = "data", assay="SCT")
  
  filtered_data <- sparsityFiltering(data, sample_thresh = 1000, feature_thresh = 700)
  
  filtered_seurat_obj <- subset(seurat_obj, cells = colnames(filtered_data), features = rownames(filtered_data))
  
  print(paste0("subsetNIData: Sparsity filtered - ", ncol(filtered_seurat_obj), " cells and ", nrow(filtered_seurat_obj), " features"))
  
  ## 3) Subset just the cell types of interest
  Idents(filtered_seurat_obj) <- filtered_seurat_obj[[celltype_metadata_key]] %>% rownames_to_column() %>% deframe()
  filtered_seurat_obj <- subset(filtered_seurat_obj, idents = celltype_metadata_val)
  
  print(paste0("subsetNIData: Subset to cell types of interest (", celltype_metadata_key , "=", paste(unlist(comp$metadata_val), collapse = "+") , ") - ", ncol(filtered_seurat_obj), " cells and ", nrow(filtered_seurat_obj), " features"))
  
  ## 4) Then subset just those belonging to the conditions of interest
  Idents(filtered_seurat_obj) <- filtered_seurat_obj$condition
  filtered_seurat_obj <- subset(filtered_seurat_obj, idents = conditions)
  
  print(paste0("subsetNIData: Subset to conditions of interest (", paste(conditions, collapse = " vs ") , ") - ", ncol(filtered_seurat_obj), " cells and ", nrow(filtered_seurat_obj), " features"))
  
  ## 5) If we have more than the number of desired features (i.e., genes) 
  ## then variance-based filter down to that number
  if (length(Features(filtered_seurat_obj)) > max_features) {
    
    # We copy here so as to not mess up the SCT counts in the original object
    seurat_obj_copy <- filtered_seurat_obj
    
    seurat_obj_copy <- SCTransform(seurat_obj_copy, 
                                   vst.flavor = "v2", 
                                   return.only.var.genes = TRUE, 
                                   ncells = length(Cells(seurat_obj_copy)), 
                                   variable.features.n = max_features,
                                   verbose = TRUE)
    
    filtered_seurat_obj <- subset(filtered_seurat_obj, features = VariableFeatures(seurat_obj_copy))
    
    print(paste0("subsetNIData: Variance filtered - ", ncol(filtered_seurat_obj), " cells and ", nrow(filtered_seurat_obj), " features"))
  }
  
  ## 6) Subset the number of samples if needed and evenly distribute them 
  ## between the two conditions and cell types
  celltype_and_conditions_table <- table(filtered_seurat_obj[[celltype_metadata_key]] %>% rownames_to_column() %>% deframe(), 
                                         filtered_seurat_obj$condition)
  
  print(celltype_and_conditions_table)
  
  celltype_and_conditions_min <- min(celltype_and_conditions_table[celltype_and_conditions_table > 0])
  
  data <- GetAssayData(filtered_seurat_obj, layer = "data", assay="SCT") %>% 
    as.matrix() %>% t() %>% as.data.frame() %>% rownames_to_column("barcode") %>% 
    left_join(tibble(barcode = names(filtered_seurat_obj$condition), y = filtered_seurat_obj$condition), 
              by = c("barcode" = "barcode")) %>% 
    mutate(y = factor(y, levels = conditions)) %>% 
    relocate(y) %>% 
    left_join(filtered_seurat_obj[[celltype_metadata_key]] %>% rownames_to_column("barcode"), 
              by = c("barcode" = "barcode")) %>% 
    relocate(any_of(celltype_metadata_key), .after = y) %>% 
    group_by(across(all_of(c("y", celltype_metadata_key)))) %>% 
    slice_sample(n = celltype_and_conditions_min)
  
  print(table(data[[celltype_metadata_key]], data$y))
  
  data <- data %>% 
    ungroup() %>% 
    select(-any_of(celltype_metadata_key)) %>% 
    column_to_rownames("barcode")
  
  print(paste0("subsetNIData: Downsampling/Balancing across conditions/celltypes - ", nrow(data), " cells and ", ncol(data) - 1, " features"))
  
  ## 7) Now remove any features remaining with ZERO variance (SLIDE will break if it gets features with zero variance/StdDev)
  non_zero_variance_features <- data %>% select(-y) %>% as.matrix() %>% matrixStats::colVars() %>% tibble::enframe() %>% filter(value > 0) %>% pull(name)
  
  data <- data %>% select(y, all_of(non_zero_variance_features))
  
  print(paste0("subsetNIData: Removing zero-variance features - ", nrow(data), " cells and ", ncol(data) - 1, " features"))
  
  return(data)
}
  
prepSLIDE <- function(data,
                      out_path = "./",
                      lambda = c(0.1, 0.5, 1.0),
                      delta =  c(0.01, 0.1), 
                      spec = 0.3, 
                      thresh_fdr = 0.2) {
  
  yaml_args <- list()
  
  yaml_args$y_factor <- TRUE
  #yaml_args$alpha_level <- 0.05
  yaml_args$thresh_fdr <- thresh_fdr
  yaml_args$delta <- delta
  yaml_args$lambda <- lambda
  yaml_args$spec <- spec
  
  yaml_args$SLIDE_iter <- 1000
  yaml_args$SLIDE_top_feats <- 50
  yaml_args$CViter <- 10
  yaml_args$sampleCV_K <- 4
  yaml_args$do_interacts <- TRUE
  
  # Create out folder if it doesn't exist
  yaml_args$out_path <- paste0(sub("\\/$", "", out_path), "/")
  dir.create(yaml_args$out_path, showWarnings = FALSE)
  
  yaml_args$eval_type <- ifelse(length(unique(data$y)) > 2, "corr", "auc")
  
  x <- data %>% select(-y) %>% rownames_to_column("barcode")
  y <- data %>% select(y) %>% mutate(y = as.numeric(y) - min(as.numeric(y))) %>% rownames_to_column("barcode") 
  
  yaml_args$y_levels <- sort(unique(y$y))
  
  yaml_args$x_path <- paste0(yaml_args$out_path, "x.csv")
  readr::write_csv(x, yaml_args$x_path, col_names = TRUE) 
  
  yaml_args$y_path <- paste0(yaml_args$out_path, "y.csv")
  readr::write_csv(y, yaml_args$y_path, col_names = TRUE) 
  
  
  ## Write YAML file
  
  yaml_path = paste0(yaml_args$out_path, "er.yaml")
  
  yaml::write_yaml(yaml_args, yaml_path)
  
  # yaml_string <- paste0("x_path: ", yaml_args$x_path, "\n", 
  #                       "y_path: ", yaml_args$y_path, "\n", 
  #                       "out_path: ", yaml_args$out_path, "\n", 
  #                       "y_factor: ", yaml_args$y_factor, "\n", 
  #                       "y_levels: [", paste(yaml_args$y_levels, collapse = ', '), "]\n", 
  #                       "eval_type: ", yaml_args$eval_type, "\n", 
  #                       "rep_cv: ", yaml_args$rep_cv, "\n", 
  #                       "alpha_level: ", yaml_args$alpha_level, "\n", 
  #                       "thresh_fdr: ", yaml_args$thresh_fdr, "\n", 
  #                       "std_cv: ", yaml_args$std_cv, "\n", 
  #                       "std_y: ", yaml_args$std_y, "\n", 
  #                       "k: ", yaml_args$k, "\n", 
  #                       "nreps: ", yaml_args$nreps, "\n", 
  #                       "permute: ", yaml_args$permute, "\n", 
  #                       "benchmark: ", yaml_args$benchmark, "\n", 
  #                       "delta: ", yaml_args$delta, "\n", 
  #                       "lambda: ", yaml_args$lambda, "\n", 
  #                       "spec: ", yaml_args$spec, "\n")
  # readr::write_file(yaml_string, yaml_path)
  
  return(0)
}

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

