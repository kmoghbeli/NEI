#Slide Helper Functions

options(future.globals.maxSize = 1000 * 1024^8)   ## 8GB

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

subsetNIData <- function(seurat_obj, 
                         celltype_metadata_key, 
                         celltype_metadata_val, 
                         conditions, 
                         exclude_mito_ribo_genes = TRUE, 
                         genes_to_keep = NULL) {
  
  max_cells_per_condition <- 10000000
  max_features <- 5000
  
  message("subsetNIData: Initial dataset - ", ncol(seurat_obj), " cells and ", nrow(seurat_obj), " features")
  
  filtered_seurat_obj <- seurat_obj
  
  ## 1) Remove mitochondrial and ribosomal genes
  if (TRUE == exclude_mito_ribo_genes) {
    non_mt_ribo_genes <- which(!grepl("^MT-|^RPS|^RPL|^mRP|^ATP", 
                                      Features(filtered_seurat_obj), 
                                      ignore.case = TRUE))
    
    filtered_seurat_obj <- subset(filtered_seurat_obj, features = Features(filtered_seurat_obj)[non_mt_ribo_genes])
  }  
  
  message("subsetNIData: Removed mito/ribo genes - ", ncol(filtered_seurat_obj), " cells and ", nrow(filtered_seurat_obj), " features")
  
  ## 2) First Zero Filter
  data <- GetAssayData(filtered_seurat_obj, layer = "data", assay="SCT")
  
  filtered_data <- sparsityFiltering(data, sample_thresh = 1000, feature_thresh = 700)
  
  filtered_seurat_obj <- subset(filtered_seurat_obj, cells = colnames(filtered_data), features = rownames(filtered_data))
  
  message("subsetNIData: Sparsity filtered - ", ncol(filtered_seurat_obj), " cells and ", nrow(filtered_seurat_obj), " features")
  
  ## 3) Subset just the cell types of interest
  Idents(filtered_seurat_obj) <- filtered_seurat_obj[[celltype_metadata_key]] %>% rownames_to_column() %>% deframe()
  filtered_seurat_obj <- subset(filtered_seurat_obj, idents = celltype_metadata_val)
  
  message("subsetNIData: Subset to cell types of interest (", celltype_metadata_key , "=", paste(unlist(celltype_metadata_val), collapse = "+") , ") - ", ncol(filtered_seurat_obj), " cells and ", nrow(filtered_seurat_obj), " features")
  
  ## 4) Then subset just those belonging to the conditions of interest
  Idents(filtered_seurat_obj) <- filtered_seurat_obj$condition
  filtered_seurat_obj <- subset(filtered_seurat_obj, idents = conditions)
  
  message("subsetNIData: Subset to conditions of interest (", paste(conditions, collapse = " vs ") , ") - ", ncol(filtered_seurat_obj), " cells and ", nrow(filtered_seurat_obj), " features")
  
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
    
    message("subsetNIData: Variance filtered - ", ncol(filtered_seurat_obj), " cells and ", nrow(filtered_seurat_obj), " features")
  }
  
  ## 6) Subset the number of samples if needed and evenly distribute them 
  ## between the two conditions and cell types
  celltype_and_conditions_table <- table(filtered_seurat_obj[[celltype_metadata_key]] %>% rownames_to_column() %>% deframe(), 
                                         filtered_seurat_obj$condition)
  
  #message(celltype_and_conditions_table)
  
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
  
  #message(table(data[[celltype_metadata_key]], data$y))
  
  data <- data %>% 
    ungroup() %>% 
    select(-any_of(celltype_metadata_key)) %>% 
    column_to_rownames("barcode")
  
  message("subsetNIData: Downsampling/Balancing across conditions/celltypes - ", nrow(data), " cells and ", ncol(data) - 1, " features")
  
  ## 7) Now remove any features remaining with ZERO variance (SLIDE will break if it gets features with zero variance/StdDev)
  non_zero_variance_features <- data %>% select(-y) %>% as.matrix() %>% matrixStats::colVars() %>% tibble::enframe() %>% filter(value > 0) %>% pull(name)
  
  data <- data %>% select(y, all_of(non_zero_variance_features))
  
  message("subsetNIData: Removing zero-variance features - ", nrow(data), " cells and ", ncol(data) - 1, " features")
  
  ## 8) Add back in any genes in the passed in "genes_to_keep" list that we may have removed in our filtering steps above
  if (!is.null(genes_to_keep) && !is.na(genes_to_keep) && length(genes_to_keep) > 0) {
    genes_to_keep <- setdiff(genes_to_keep, colnames(data))   # Only add the ones we don't already have
    
    counts_for_genes_to_keep <- 
      GetAssayData(seurat_obj, layer = "data", assay="SCT")[genes_to_keep, rownames(data)] %>% 
      as.matrix() %>% t() %>% as.data.frame() %>% rownames_to_column("barcode")
    
    data <- left_join(x = data %>% rownames_to_column("barcode"), 
                      y = counts_for_genes_to_keep, 
                      by = join_by("barcode" == "barcode")) %>% 
      column_to_rownames("barcode")
  }
  
  
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
  
  yaml_path = paste0(yaml_args$out_path, "slide.yaml")
  
  yaml::write_yaml(yaml_args, yaml_path)
}