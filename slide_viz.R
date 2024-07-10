## Visualization Functions for SLIDE results and the identified Latent Factors

library(tidyverse)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(gprofiler2)
library(qgraph)
conflicted::conflict_prefer_all("dplyr", quiet = TRUE)


kaveh_colors <- c("#40004b", "#c2a5cf", "#f5f5f5", "#80cdc1", "#018571")

## Covariance Graphs
plot_lf_covar_graphs <- function(config_yaml) {
  for (result in config_yaml$er_results) {
    
    er_x <- readr::read_csv(paste0(result$path, "x.csv"), show_col_types = FALSE)
    
    lfs <- str_split_1(result$lfs, ",[ ]?")
    
    for (lf in lfs) {
      
      lf_features <- readr::read_table(paste0(result$path, "gene_list_Z", lf, ".txt"), show_col_types = FALSE) %>% 
        mutate(names = ifelse(grepl("^X.+Rik$", names), sub("^X", "", names), names)) %>% 
        mutate(names = gsub("\\.", "-", names))   # For some reason ER/SLIDE results swap out "-" for "." in our gene names
      
      corr_matrix <- cor(er_x %>% select(all_of(lf_features %>% pull(names))))
      
      # In ER/SLIDE results: Red genes are high in 1/Experimental, blue genes are high in 0/Control
      
      ## [BEG] Covar Graph ##
      pdf(file = paste0(result$path, "Z", lf, "_covar.pdf"))
      
      qgraph(
        corr_matrix,
        layout = "spring",
        #threshold = 0.5,
        threshold = 0,
        repulsion = 0.1, 
        node.width = 1.4, 
        node.height = 1.4, 
        labels = colnames(corr_matrix),
        label.font = 4,
        label.cex = 2,
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
      ## [END] Covar Graph ##
    }
  }
}


## LF vs LF Scatter Plots
plot_lf_scatterplots <- function(config_yaml) {
  for (result in config_yaml$er_results) {
    
    y <- readr::read_csv(paste0(result$path, "y.csv"), show_col_types = FALSE) %>% 
      mutate(condition = ifelse(0 == y, result$condition1, result$condition2))
    
    z_mat <- readr::read_csv(paste0(result$path, "z_matrix.csv"), show_col_types = FALSE) %>% 
      rename(barcode = `...1`) %>% 
      left_join(y, by = join_by(barcode == barcode))
    
    
    lfs <- paste0("Z", str_split_1(result$lfs, ",[ ]?"))
    
    lf_pairs <- combn(lfs, 2) %>% t()
    
    for (i in 1:nrow(lf_pairs)) {
      
      p <- z_mat %>% 
        ggplot(aes(.data[[lf_pairs[i, 1]]], .data[[lf_pairs[i, 2]]], color = condition)) + 
        ggprism::scale_color_prism("colorblind_safe") + 
        geom_point() + 
        ggpubr::theme_pubclean() + 
        ggtitle(paste0(result$title, " (", result$condition1, " vs ", result$condition2, ")"))
      
      ggsave(file = paste0(result$path, lf_pairs[i, 1], "_vs_", lf_pairs[i, 2], "_scatter.pdf"), height = 7, width = 7)
    }
  }
}


## Pathways Analysis
### Example custom GMT uploaded to GProfiler (e.g., to use MSigDB GMT files)
### Can get them here: https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp
## Canonical Pathways sets
## upload_GMT_file(gmtfile = "c2.cp.v2023.1.Hs.symbols.gmt")
# Your custom annotations ID is gp__deB5_Np8D_xYo
# You can use this ID as an 'organism' name in all the related enrichment tests against this custom source.
# Just use: gost(my_genes, organism = 'gp__deB5_Np8D_xYo')

plot_pathways_analysis <- function(config_yaml) {
  
  # msig_IDs <- data.frame(id = c("gp__NwdD_iLvQ_QbY",   # Hallmark Mouse
  #                               "gp__hIhS_hLOK_V4Q"),     # Canonical Pathways Mouse, 
  #                        name = c("MSigDB Hallmark",
  #                                 #"MSigDB BioCarta", 
  #                                 #"MSigDB ImmuneSigDB", 
  #                                 "MSigDB Canonical"))
  
  for (result in config_yaml$er_results) {
    
    lfs <- str_split_1(result$lfs, ",[ ]?")
    
    pooled_lf_genes <- c()
    
    for (lf in lfs) {
      
      lf_features <- readr::read_table(paste0(result$path, "gene_list_Z", lf, ".txt"), show_col_types = FALSE) %>% 
        mutate(names = ifelse(grepl("^X.+Rik$", names), sub("^X", "", names), names)) %>% 
        mutate(names = gsub("\\.", "-", names))   # For some reason ER/SLIDE results swap out "-" for "." in our gene names
      
      pooled_lf_genes <- c(pooled_lf_genes, lf_features$names)
      
      filename <- paste0("./", result$path, "Z", lf, "_pathways.pdf")
      
      gprofiler_result <- gprofiler2::gost(lf_features$names, 
                                           organism="mmusculus", 
                                           ordered_query = FALSE,  # True for GSEA style p-values, False for standard ORA
                                           user_threshold = 0.05, 
                                           correction_method = "fdr")
      
      p <- pathways_plot(gprofiler_result, paste0("Z", lf, " ORA"))
      
      ggsave(filename, plot = p, width = 12, height = 9, units = "in") 
    }
    
    # Run the Pooled ORA for the SLIDE result LFs
    gprofiler_result <- gprofiler2::gost(unique(pooled_lf_genes), 
                                         organism="mmusculus", 
                                         ordered_query = FALSE,  # True for GSEA style p-values, False for standard ORA
                                         user_threshold = 0.05, 
                                         correction_method = "fdr")
    
    p <- pathways_plot(gprofiler_result, "Pooled ORA")
    
    ggsave(paste0("./", result$path, "pooled_pathways.pdf"), 
           plot = p, width = 12, height = 9, units = "in") 
    
  }
}

pathways_plot <- function(gprofiler_result, plot_title = "") {
  
  plot <- 
  gprofiler_result$result %>%
    filter(grepl("GO:BP", source) | grepl("KEGG", source) | grepl("REAC", source)) %>%
    mutate(term_name = stringr::str_trunc(term_name, 50),
           neg_log_p = -log(p_value), 
           percent_intersect = (intersection_size / term_size) * 100) %>% 
    dplyr::select(term_id, source, term_name, neg_log_p, percent_intersect) %>% 
    group_by(source) %>% 
    arrange(desc(source), desc(neg_log_p), desc(percent_intersect)) %>% 
    slice_head(n = 10) %>% 
    ggplot(aes(x = neg_log_p, y = reorder(term_name, neg_log_p, decreased = TRUE), fill = source)) + 
    geom_bar(stat = "identity") + 
    ggpubr::theme_pubclean() + 
    scale_fill_viridis_d() + 
    labs(x = "-log10p", y = element_blank(), title = plot_title) + 
    theme(plot.title = element_text(hjust = 0.5, face = "bold"), 
          axis.text = element_text(face = "bold", size = 18),
          legend.position = "bottom", 
          legend.direction = "horizontal")
  
  return(plot)
}


## Below is the expected format for sample Config YAML ##
config <- yaml::yaml.load("
er_results:

  # - result:
  #   path: 'slide_runs/combined_tg_Mac_KOSvRE_D0.01_L1_S0.5/'
  #   lfs: '18, 27, 31, 107, 147, 157'
  #   condition1: 'KOS'
  #   condition2: 'RE'
  #   title: 'Combined TG Macs'

  - result:
    path: 'slide_runs/combined_tg_Mac_KOSvRE_D0.01_L0.5_S0.3/'
    lfs: '87, 27, 139, 18, 147, 31, 107, 157'
    condition1: 'KOS'
    condition2: 'RE'
    title: 'Combined TG Macs'
")


plot_lf_covar_graphs(config)
plot_lf_scatterplots(config)
plot_pathways_analysis(config)

## Transcription Factor stuff
# stat1_targets <- jsonlite::fromJSON("stat1_targets.json")
# stat1_targets$associations$gene$symbol
# 
# 
# ## SLIDE results exploration
# x <- readr::read_csv("slide_runs/combined_tg_T_KOSvRE_D0.01_L1_S0.5/x.csv")
# y <- readr::read_csv("slide_runs/combined_tg_T_KOSvRE_D0.01_L1_S0.5/y.csv")
# 
# controlPerf <- readr::read_rds("slide_runs/OLD_combined_tg_Mac_KOSvRE_D0.01_L1.0_S0.5/ControlPerformance.rds")
# controlPerfPlot <- readr::read_rds("slide_runs/OLD_combined_tg_Mac_KOSvRE_D0.01_L1.0_S0.5/ControlPerformancePlot.rds")
# 
# final_delta <- readr::read_rds("slide_runs/combined_tg_T_KOSvRE_D0.01_L1_S0.5/final_delta_0.01_lambda_0.1.rds")
# plotSigGenes <- readr::read_rds("slide_runs/combined_tg_NK_KOSvRE_D0.01_L0.1/plotSigGenes_data.RDS")
# pipeline3_y_mapping <- readr::read_rds("slide_runs/combined_tg_NK_KOSvRE_D0.01_L0.1/pipeline3_y_mapping.rds")


# base_slide_path <- "slide_runs/combined_tg_Mac_KOSvRE_D0.01_L0.5_S0.3/"
# lfs <- c(87, 27, 139, 18, 147, 31, 107, 157)
# lf_features <- tibble()
# 
# for (lf in lfs) {
#   
#   lf_features <- 
#     bind_rows(lf_features, 
#               readr::read_table(paste0(base_slide_path, "gene_list_Z", lf, ".txt"), show_col_types = FALSE) %>% 
#                 mutate(names = ifelse(grepl("^X.+Rik$", names), sub("^X", "", names), names)) %>% 
#                 mutate(names = gsub("\\.", "-", names)) %>%  # For some reason ER/SLIDE results swap out "-" for "." in our gene names 
#                 mutate(lf_num = as.numeric(lf)) %>% 
#                 rename(gene = names)     
#     )
# }
# 
# plotSigGenes <- readr::read_rds("slide_runs/combined_tg_Mac_KOSvRE_D0.01_L0.5_S0.3/plotSigGenes_data.RDS") %>% 
#   mutate(names = ifelse(grepl("^X.+Rik$", names), sub("^X", "", names), names)) %>% 
#   mutate(names = gsub("\\.", "-", names)) %>%  # For some reason ER/SLIDE results swap out "-" for "." in our gene names
#   filter(lf_num %in% lfs)
# 
# blah <- lf_features %>% inner_join(plotSigGenes, by = join_by(gene == names, lf_num == lf_num)) %>% 
#   select(gene, A_loading.x, A_loading.y, AUCs.x, AUCs.y, corrs.x, corrs.y, color.x, color.y)
