## Covariance Graphs

library(tidyverse)
library(EssReg)
library(SLIDE)
library(SLIDEHelper)
library(qgraph)

kaveh_colors <- c("#40004b", "#c2a5cf", "#f5f5f5", "#80cdc1", "#018571")

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

config <- yaml::yaml.load("
er_results:
  - result:
    path: 'july_cornea_T_CONTROLvRE_D0.01_L0.1/'
    lfs: '3, 17, 20'
    condition1: 'Control'
    condition2: 'RE'
    title: 'July Cornea T Cells'

  - result:
    path: 'july_cornea_Mye_CONTROLvRE_D0.01_L0.1/'
    lfs: '19, 29, 31'
    condition1: 'Control'
    condition2: 'RE'
    title: 'July Cornea Myeloid Cells'

  - result:
    path: 'july_cornea_Epi_CONTROLvRE_D0.01_L0.1/'
    lfs: '5, 34, 38'
    condition1: 'Control'
    condition2: 'RE'
    title: 'July Cornea Epithelial Cells'
")

for (result in config$er_results) {

  er_x <- readr::read_csv(paste0(result$path, "x.csv"), show_col_types = FALSE)
  
  lfs <- str_split_1(result$lfs, ",[ ]?")
  
  for (lf in lfs) {
    
    lf_features <- readr::read_table(paste0(result$path, "gene_list_Z", lf, ".txt"), show_col_types = FALSE) %>% 
      mutate(names = ifelse(grepl("^X.+Rik$", names), sub("^X", "", names), names)) %>% 
      mutate(names = gsub("\\.", "-", names))   # For some reason ER/SLIDE results swap out "-" for "." in our gene names
    
    corr_matrix <- cor(er_x %>% select(all_of(lf_features %>% pull(names))))
    
    # In ER/SLIDE results: Red genes are high in 1/Experimental, blue genes are high in 0/Control
    
    pdf(file = paste0(result$path, "Z", lf, "_covar.pdf"))
    
    qgraph(
      corr_matrix,
      layout = "spring",
      #threshold = 0.5,
      threshold = 0,
      repulsion = 0.1,
      labels = colnames(corr_matrix),
      label.font = 4,
      label.scale.equal = TRUE,
      label.prop = 0.95,
      shape = ifelse("Blue" == lf_features$color, "ellipse", "rectangle"),  
      color = ifelse("Blue" == lf_features$color, "#8EE5EE", "#FFB5C5"),
      posCol = kaveh_colors[length(kaveh_colors)],
      negCol = kaveh_colors[1],
      #height = 5,
      #width = 7,
      label.cex = 1.25  
    )
    
    legend("topright", legend = c("Positive Correlation", "Negative Correlation"), 
           title = "Edge Colors", col = c(kaveh_colors[length(kaveh_colors)], kaveh_colors[1]), lty = 1, lwd = 3.5, cex = 0.7)
    legend("topleft", legend = c(paste0("High in ", result$condition1), paste0("High in ", result$condition2)), 
           title = "Genes", pch = c(19, 15), col = c("#8EE5EE", "#FFB5C5"), pt.cex = 2, cex = 0.7)
    title(paste0("Z", lf), line = 3)
    
    dev.off()
  }
}










# 
# 
# 
# 
# ######## [END] ZARIFEH CODE #########
# 
# ######## [BEG] MARISA CODE #########
# library(ggplot2) 
# library(dplyr)
# library(tidyr)
# library(stringr)
# library(qgraph) 
# # Define a named vector that associates cell types with colors
# cellTypeColors <- c(
#   "BC" = "#888888",
#   "PB" = "#E69F00",
#   "CD8T" = "#56B4E9",
#   "CD4Tconv" = "#009E73",
#   "CD4Treg" = "#40E0D0",
#   "CD14CD16macro" = "#D55E00",
#   "CD14mono" ="#0072B2", 
#   "CD16mono" = "#D2B48C"  ,   
#   "CD14CD16mono" =  "#D7BFDC",  
#   "pDCs" = "#F0E442",         
#   "NK" =  "#CC79A7"       
# )
# ############################################################ 
# # plotShape function
# ############################################################ 
# plotShape <- function(yaml_path, lf, threshold, repulsion) {
#   input <- yaml::yaml.load_file(yaml_path)
#   x <- read.csv(input$x_path, row.names = 1, header = TRUE)
#   
#   sigGenes <- readRDS(paste0(input$out_path, "plotSigGenes_data.RDS"))
#   
#   # Get specific lf 
#   lf_gene_df <- sigGenes[sigGenes$lf_num == lf,]
#   
#   # Unique df
#   lf_gene_df <-  lf_gene_df %>%
#     distinct(names, .keep_all = TRUE)
#   
#   # Unique lf gene names & color
#   lf_genes <- lf_gene_df$names
#   lf_gene_colors <- lf_gene_df$color
#   lf_shape = recode(lf_gene_colors, 'Blue'='circle', .default = 'square')
#   
#   # correlation matrix for lf 
#   sub <- x[colnames(x) %in% lf_genes]
#   corr_matrix <- cor(sub)
#   
#   # cell type and gene info 
#   cell_type <- str_extract(colnames(corr_matrix), "(?<=_).*")
#   group <- factor(cell_type)
#   gene <-  str_extract(colnames(corr_matrix), ".*(?=_)")
#   
#   q <- qgraph(
#     corr_matrix,
#     layout = "spring",  
#     #groups = group,
#     colors = cellTypeColors[cell_type],  # Use the defined colors based on cell type
#     shape = lf_shape, 
#     rescale = T,
#     labels = gene,
#     border.width = 2.5, 
#     posCol="#B0C4DE",
#     negCol=  "#C70039",
#     fade = T,
#     alpha = 0.05,
#     threshold = threshold,
#     repulsion = repulsion,
#     label.prop = 0.9
#   ) 
#   q
# }
# ### EXAMPLE ############ 
# # plots latent factor 2, threshold = 0.3, and repulsion = 0.75
# #cr_paths= c( 'best_models/CR/pbmc/cr_0.05_1.yaml', 'best_models/CR/til5ct/cr_0.05_0.5.yaml', 'best_models/CR/til_allct/cr_0.1_1.yaml')
# #plotShape(cr_paths[2], 2, 0.3, 0.75) 
# 
# ######## [END] MARISA CODE #########