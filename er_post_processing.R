library(tidyverse)
library(doParallel)

plot_slide_violin <-function(lambda_rep) {
  p <- 
  lambda_rep %>% 
    filter(method %in% c("plainER", "plainER_y")) %>% 
    mutate(method = sub("plain", "", method), 
           method = sub("_y", "\n(perm)", method), 
           method = sub("ER", "SLIDE", method),
           method = factor(method)) %>% 
    rename(any_of(c(res = "auc", res = "corr"))) %>% 
    ggplot(aes(x = method,
               y = res,
               fill = method)) +
    ggplot2::geom_violin() +
    ggpubr::theme_pubr() + 
    labs(y = "AUC") + 
    ggprism::scale_fill_prism(palette = "colorblind_safe") + 
    theme(axis.title.x = element_blank(), 
          axis.text = element_text(size = 20, face = "bold"), 
          axis.text.x = element_text(angle = 0, hjust = 0.5), 
          title = element_text(size = 24, face = "bold"), 
          axis.line =  element_line(colour = 'black', linewidth = 2), 
          #axis.ticks = element_line(colour = 'black', linewidth = 2), 
          legend.position = "none", 
          aspect.ratio = 1.7)
  
  return(p)
}

registerDoParallel(detectCores() - 2)

########################################################################

config <- yaml::yaml.load("
er_results: 
  #- path: 'slide_runs/combined_tg_Mac_KOSvRE_D0.01_L1_S0.5/' 
  - path: 'slide_runs/combined_tg_Mac_KOSvRE_D0.01_L0.5_S0.3/'
")

for (result in config$er_results) {
  
  print(paste0("Processing: ", result$path))
  
  yaml_path <- paste0(result$path, "er.yaml")
  
  yaml_args <- yaml::yaml.load_file(yaml_path)
  
  # Set a default "spec"
  yaml_args$spec <- ifelse(is.null(yaml_args$spec), 0.1, yaml_args$spec)
  
  # check for ER results in output folder
  er_results_path = list.files(yaml_args$out_path, pattern = "final_delta")
  
  if (0 == length(er_results_path)) {
    print("ER did not complete successfully.")
    next
  }
  
  lambda_df <- readr::read_rds(paste0(yaml_args$out_path, "pipeline_step5.rds"))
  
  er_results_plot <- plot_slide_violin(lambda_df)
  
  ggsave("slide_auc_corr.png", width = 4, height = 5, path = yaml_args$out_path)
  
  
  er_results_path <- paste0(yaml_args$out_path, er_results_path)
  
  Z_matrix <- SLIDEHelper::CalcZMatrix(yaml_args$x_path, er_results_path, yaml_args$out_path)
  
  # run slide on these results
  SLIDE_res = SLIDEHelper::runSLIDE(y_path = yaml_args$y_path,
                                    z_path = NULL,
                                    z_matrix = Z_matrix,
                                    er_path = er_results_path,
                                    do_interacts = TRUE,
                                    spec = yaml_args$spec,
                                    niter = 100)
  
  n_feats <- 10
  condition <- yaml_args$eval_type
  
  SLIDE_res <- SLIDEHelper::GetTopFeatures(yaml_args$x_path, 
                                           yaml_args$y_path, 
                                           er_results_path, 
                                           yaml_args$out_path,
                                           SLIDE_res, 
                                           num_top_feats = n_feats, 
                                           condition)
  
  sig_gene_plots <- SLIDE::plotSigGenes(SLIDE_res, plot_interaction = TRUE, yaml_args$out_path)
  
  sig_gene_plots %>% readr::write_rds(paste0(yaml_args$out_path, "sig_gene_plots.rds"))
  
  y_matrix <- readr::read_csv(yaml_args$y_path) %>% column_to_rownames("barcode")
  
  SLIDE::calcControlPerformance(z_matrix = Z_matrix,
                                y = y_matrix, 
                                do_interacts = TRUE,
                                SLIDE_res = SLIDE_res,     # the SLIDE_res has to be the output from GetTopFeatures
                                condition = condition,
                                out_path = yaml_args$out_path)
  
  controlPerfPlot <- readr::read_rds(paste0(yaml_args$out_path, "ControlPerformancePlot.rds"))
  ggsave(filename = "ControlPerformancePlot.png", plot = controlPerfPlot, path = yaml_args$out_path, width = 7, height = 5)
}
