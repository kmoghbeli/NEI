# Seurat Pipeline

library(tidyverse)
library(ggpubr)
library(ggprism)
library(Seurat)
library(harmony)
library(leiden)
library(scRepertoire)
conflicted::conflict_prefer_all("dplyr", quiet = TRUE)

options(Seurat.object.assay.version = "v3")

#reticulate::py_install("leidenalg")
#leidenalg <- reticulate::import("leidenalg")


# Parallell Processing
#future::plan(future::multisession, workers = future::availableCores() - 2)
num_cores <- future::availableCores()
if(num_cores > 3) {
  library(doParallel)
  
  cl <- makePSOCKcluster(num_cores - 2)
  registerDoParallel(cl)
}


seurat_norm_integrate <- function(config, seurat_verbose = FALSE) {
 
  # Set reasonable defaults
  config$normalization <- ifelse(is.null(config$normalization), "sct", config$normalization)
  config$seurat_assay <- ifelse("sct" == config$normalization, "SCT", "RNA")
  
  seurat_objs <- list()
  
  cat("Processing: ", config$filename, "\n")
  
  for (sample in config$samples) {
    
    seurat_obj <- NULL
    
    sample_counts <- Read10X_h5(paste0(sample$path, "/filtered_feature_bc_matrix.h5"))
    
    if (is.null(sample$htos)) {
      ### No HTOs
      seurat_obj <- CreateSeuratObject(sample_counts, project = config$project, min.cells = 3, min.features = 200)
      
    } else {
      
      ### Yes HTOs
      
      ## Step 1: Prep the HTO data for future de-multiplexing
      # Select cell barcodes detected by both RNA and HTO In the example datasets we have already
      # filtered the cells for you, but perform this step for clarity.
      joint_barcodes <- intersect(colnames(sample_counts$`Gene Expression`),   # UMI barcodes
                                  colnames(sample_counts$`Antibody Capture`))  # HTO barcodes
      
      umis <- sample_counts$`Gene Expression`[, joint_barcodes]
      
      seurat_obj <- CreateSeuratObject(counts = umis, project = config$project, min.cells = 3, min.features = 200)
      
      htos <- sample_counts$`Antibody Capture`[, Cells(seurat_obj)]
      
      # Add HTO data as a new assay independent from RNA
      seurat_obj[["HTO"]] <- CreateAssayObject(counts = htos)
      
      # Normalize HTO data, here we use centered log-ratio (CLR) transformation
      # (We will use this later for demultiplexing)
      seurat_obj <- NormalizeData(seurat_obj, assay = "HTO", normalization.method = "CLR")
    }
    
    
    ## Remove any genes we don't want to process here
    ## (Probably not a very efficient method)
    if (!is.null(config$remove_features)) {
      seurat_obj <- subset(seurat_obj, features = setdiff(Features(seurat_obj), config$remove_features))
    }
    
    
    ## Load in and attach any TCR annotations
    if (!is.null(sample$tcr_path)) {
      sample_contigs <- readr::read_csv(paste0(sample$tcr_path, "/filtered_contig_annotations.csv"))
      
      #print(sample_contigs %>% group_by(barcode) %>% summarise(count = n()), n = Inf)  # just for debugging
      
      tcr_data <- combineTCR(sample_contigs, removeNA = FALSE, removeMulti = FALSE, filterMulti = TRUE)
      
      seurat_obj <- combineExpression(tcr_data,
                                      seurat_obj,
                                      cloneCall = "aa", 
                                      chain = "both",
                                      proportion = FALSE, 
                                      cloneSize = c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))
    }
    
    
    # Initial QC (mitochondrial genes, low features, doublets)
    seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern="^MT-|^mt-")
    
    ## Sample-specific thresholding/filtering: https://github.com/satijalab/seurat/issues/3396
    print(VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
    
    feature_max <- quantile(seurat_obj$nFeature_RNA, probs = 0.999)
    count_max <- quantile(seurat_obj$nCount_RNA, probs = 0.999)
    
    # percent_mito <- ifelse(is.null(config$percent_mito), 
    #                        quantile(seurat_obj$percent.mt, probs = 0.99), 
    #                        config$percent_mito)
    
    percent_mito <- ifelse(is.null(config$percent_mito), 
                           10,
                           config$percent_mito)
    
    cat("Using QC filter cutoffs: nFeature_RNA > 200 & nFeature_RNA <", feature_max, "& nCount_RNA <", count_max, "& Percent mito <", percent_mito ,"\n")
    cat("Before QC filtering:", nrow(seurat_obj), "genes x", ncol(seurat_obj), "cells\n\n")
    
    seurat_obj <- subset(seurat_obj, 
                         subset = nFeature_RNA > 200 & 
                           nFeature_RNA < feature_max & 
                           nCount_RNA < count_max & 
                           percent.mt < percent_mito)
    
    cat("After QC filtering:", nrow(seurat_obj), "genes x", ncol(seurat_obj), "cells\n\n")

    ## Add Sample-level metadata
    for (key in names(sample)) {
      if (!(key %in% c("path", "tcr_path", "htos")) & !is.null(sample[[key]])) { 
        seurat_obj <- AddMetaData(seurat_obj, sample[[key]], key)
      }
    }
    
    
    ## SCT Normalize
    seurat_obj <- SCTransform(seurat_obj, 
                              vst.flavor = "v2",
                              vars.to.regress = "percent.mt",
                              return.only.var.genes = FALSE,
                              verbose = seurat_verbose)
    
    
    if (is.null(sample$htos)) {
      
      # Add to the list of seurat_objs
      seurat_objs[[sample$id]] <- seurat_obj
      
    } else {
      
      ##########################################
      #### Demultiplexing ######################
    
      # If you have a very large dataset we suggest using k_function = 'clara'. This is a k-medoid
      # clustering function for large applications You can also play with additional parameters (see
      # documentation for HTODemux()) to adjust the threshold for classification Here we are using
      # the default settings
      seurat_obj <- HTODemux(seurat_obj, assay = "HTO", positive.quantile = 0.99)
      
      #print(table(seurat_obj$HTO_classification.global, seurat_obj$HTO_maxID))
      
      ##*** Subset the Seurat Objects into their respective singlet HTOs
      Idents(seurat_obj) <- "HTO_classification.global"
      seurat_obj.singlets <- subset(seurat_obj, idents = "Singlet")
      Idents(seurat_obj.singlets) <- "HTO_classification"
      
      for (hto_config in sample$htos) {
        
        seurat_obj.singlets.hto <- subset(seurat_obj.singlets, idents = hto_config$hto)
        
        ## Add any HTO-specific metadata here
        for (key in names(hto_config)) {
          if (!(key %in% c("hto")) & !is.null(hto_config[[key]])) { 
            seurat_obj.singlets.hto <- AddMetaData(seurat_obj.singlets.hto, hto_config[[key]], key)
          }
        }
        
        seurat_objs[[paste0(sample$id, "_", hto_config$hto)]] <- seurat_obj.singlets.hto
      }
      
    }
  }
  
  
  integ.features <- SelectIntegrationFeatures(object.list = seurat_objs, nfeatures = 3000)
  
  ## Merge
  merged_seurat <- merge(x = seurat_objs[[1]],
                         y = seurat_objs[2:length(seurat_objs)],
                         merge.data = TRUE)
  
  VariableFeatures(merged_seurat) <- integ.features
  
  merged_seurat <- RunPCA(object = merged_seurat, assay = "SCT", npcs = 50, verbose = seurat_verbose)
  
  harmony_vars <- ifelse(is.null(config$harmony_vars), "date", config$harmony_vars)
  
  merged_seurat <- RunHarmony(object = merged_seurat,
                              assay.use = "SCT",
                              reduction = "pca",
                              dims.use = 1:50, 
                              lambda = NULL,  # Ridge regression penalty - when set to NULL, harmony tries to estimate
                              kmeans_init_nstart=20, 
                              kmeans_init_iter_max=100, 
                              group.by.vars = harmony_vars, 
                              reduction.save = "harmony", 
                              plot_convergence = seurat_verbose)
  
  #res <- try({FindClusters(object = merged_seurat, resolution = 1.0, method = "igraph", algorithm = "Leiden")})
  #res <- try({FindClusters(object = merged_seurat, resolution = 1.0, algorithm = "Leiden")})
  
  merged_seurat <- RunUMAP(object = merged_seurat, assay = "SCT", reduction = "harmony", dims = 1:50)
  merged_seurat <- FindNeighbors(object = merged_seurat, assay = "SCT", reduction = "harmony", dims = 1:50, verbose = seurat_verbose)
  
  # Cluster using Leiden algorithm if it's installed and working, otherwise use the default Louvain
  merged_seurat <- tryCatch(
    {FindClusters(object = merged_seurat, resolution = 1.0, method = "igraph", algorithm = "Leiden")}, 
    error = function(err) {
      message("Unable to run Leiden clustering, defaulting to Louvain.\n", err) 
      FindClusters(object = merged_seurat, resolution = 1.0)
    }
  )
  
  merged_seurat <- PrepSCTFindMarkers(object = merged_seurat, verbose = seurat_verbose)
  
  #merged_seurat[["RNA"]] <- as(object = merged_seurat[["RNA"]], Class = "Assay")
  
  merged_seurat %>% SeuratDisk::SaveH5Seurat(paste0(data_dir, config$filename, ".h5Seurat"), overwrite = TRUE, verbose = seurat_verbose)
}
