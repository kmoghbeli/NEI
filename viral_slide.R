# NeuroImmune SLIDE

# Initialization ----
#devtools::install_github("jishnu-lab/SLIDE")

library(tidyverse)
library(Seurat)
library(SLIDE)
conflicted::conflict_prefer_all("dplyr", quiet = TRUE)

source("slide_helper_functions.R")

data_dir <- "../data_objects/"
figures_dir <- "../figures/"


# Setup Parallel Processing ----
num_cores <- as.numeric(Sys.getenv('SLURM_CPUS_PER_TASK', unset=NA))

if(is.na(num_cores)) num_cores <- parallel::detectCores()

if(num_cores > 3) {
  library(doParallel)
  
  num_cores <- num_cores - 2
  
  registerDoParallel(cores = num_cores)
  message('Using ', num_cores, ' cores...\n')
}


dataset <- "combined_tg_viral"
celltype <- "Mye"
conditions <- c("kos", "re")
include_viral_genes <- FALSE
optional_info <- ifelse(include_viral_genes, "", "no_viral_genes")

id <- str_to_lower(paste0(dataset, "_", optional_info, "_", celltype, "_", conditions[1], "_vs_", conditions[2]))

slide_data_path = paste0("./slide_runs/", id, "/")
run_optimize_slide = !file.exists(paste0(slide_data_path, "/summary_table.csv"))

#slide_data_path <- commandArgs(trailingOnly=TRUE)

if (run_optimize_slide) {
  # Prepare Seurat data and YAML file ----
  message("Loading Seurat data [", dataset, "]...")
  
  seurat_obj <- SeuratDisk::LoadH5Seurat(paste0(data_dir, "seurat/", dataset, "_annot.h5Seurat"))
  
  hsv_genes <- rownames(seurat_obj)[str_starts(rownames(seurat_obj), "HSV-")]
  
  if (include_viral_genes) {
    genes_to_keep <- hsv_genes
  } else {
    genes_to_keep <- NULL
  }
  
  message("Subsetting Seurat data...")
  data_subset <- subsetNIData(seurat_obj, 
                              celltype_metadata_key =  "cell_L2", 
                              celltype_metadata_val = celltype, 
                              conditions = c("kos", "re"), 
                              exclude_mito_ribo_genes = TRUE,
                              genes_to_keep)
  
  message("Preparing SLIDE folder [", slide_data_path, "]...")
  prepSLIDE(data = data_subset, 
            out_path = slide_data_path,
            lambda = c(0.1, 0.5, 1.0),   # higher = sparser LFs (i.e., more loadings are zero) (0.5 - 1.0)
            delta =  c(0.01, 0.1),  # higher = fewer LFs returned from the unsupervised part (0.01, 0.1)
            spec = 0.3,     # higher = fewer significant LFs returned
            thresh_fdr = 0.2)
  

  # Run optimizeSLIDE ----
  input_params <- yaml::read_yaml(paste0(slide_data_path, "slide.yaml"))
  
  message("Checking YAML file [", slide_data_path, "slide.yaml", "]...")
  SLIDE::checkDataParams(input_params)
  
  message("Running optimizeSLIDE...")
  
  optimizeSLIDE(input_params, sink_file = FALSE)
}
  
# Check optimizeSLIDE results and run SLIDEcv to get exact performance of SLIDE model ----
slide_model_found <- FALSE
optimize_slide_results <- NULL

if (file.exists(paste0(slide_data_path, "/summary_table.csv"))) {
  optimize_slide_results <- read.csv(paste0(slide_data_path, "/summary_table.csv"), row.names = 1)

  slide_model_found = !all(is.na(optimize_slide_results$Num_of_Sig_LFs))
}

final_slide_result_exists <- FALSE
best_model <- NULL

if (slide_model_found) {
  
  # Get params from best model
  best_model <-  optimize_slide_results %>% 
    arrange(desc(sampleCV_Performance)) %>% 
    slice_head(n = 1)
  
  final_slide_result_exists <- file.exists(paste0(slide_data_path, best_model$delta, "_", best_model$lambda, "_out/slidecv.complete"))
  
} else  {
  message("No optimized slide model found:")
  message(paste0(capture.output(summary_table), collapse = "\n"))
}

if (slide_model_found && !final_slide_result_exists) {
  # create our new YAML with the optimal delta and lambda
  optimal_slide_params <- yaml::read_yaml(paste0(slide_data_path, "slide.yaml"))
  optimal_slide_params$delta <- best_model$delta
  optimal_slide_params$lambda <- best_model$lambda
  
  yaml::write_yaml(optimal_slide_params, paste0(slide_data_path, "optimal_slide.yaml"))
  
  # Run SLIDEcv
  message("Running SLIDEcv with optimized model:")
  message(paste0(capture.output(best_model), collapse = "\n"))
  
  SLIDE::SLIDEcv(paste0(slide_data_path, "optimal_slide.yaml"), 
                 nrep = 20, k = 5)
  
  # Marker to let us know that SLIDEcv was run on this param combination
  write.table(data.frame(), file = paste0(slide_data_path, best_model$delta, "_", best_model$lambda, "_out/slidecv.complete"), col.names=FALSE)
} 

if (file.exists(paste0(slide_data_path, best_model$delta, "_", best_model$lambda, "_out/slidecv.complete"))) {
  message("Final SLIDE result already exists.")
  
  optimal_slide_params <- yaml::read_yaml(paste0(slide_data_path, "optimal_slide.yaml"))
  
  SLIDE::plotCorrelationNetworks(optimal_slide_params)
}


# ## Clean up parallel stuff
# registerDoSEQ()
# 
# env <- foreach:::.foreachGlobals
# rm(list=ls(name=env), pos=env)