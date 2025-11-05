# Lasso Nested CV for Date

library(tidyverse)
library(ggpubr)
library(ggprism)
library(Seurat)
conflicted::conflict_prefer_all("dplyr", quiet = TRUE)

# DO NOT USE SEURAT V5 -> V5 objects do not play nicely yet with CellOracle in Python
options(Seurat.object.assay.version = "v3")


kaveh_colors1 <- c("#a6611a", "#dfc27d", "#f5f5f5", "#80cdc1", "#018571")
kaveh_colors2 <- c("#40004b", "#762a83", "#9970ab", "#c2a5cf", "#e7d4e8", "#f7f7f7", "#d9f0d3", "#a6dba0", "#5aae61", "#1b7837", "#00441b")
kaveh_colors3 <- c("#ff0000", "#f5f5f5","#0000ff")


config_dir <- "./config/"
data_dir <- "./data/"
model_dir <- "./models/"
results_dir <- "./results/"
figures_dir <- "./figures/"



##* Below code is for a Lasso NestedCV model
##* It expects the outcome column to be called "y"


source("ImmunoToolbox/lasso_nested_cv.R")


### LOAD DATA AND RUN THE MODEL

# Downsample
data <- subset(ni.sc,
                     subset = "Mac1" == cell_L3 &
                       ("20230629" == date | "20230824" == date) &
                       condition %in% c("control", "scratch", "kos", "re") &
                       nFeature_RNA > 500,
                     features = VariableFeatures(ni.sc),
                     downsample = 400)


### Nested CV
###------------------------------------------------------------
data_folds <- nested_cv(data, 
                        outside = vfold_cv(v = 5, repeats = 5, strata = y), 
                        inside = vfold_cv(v = 5, strata = y))

outer_results <- run_lasso_nestedcv(data_folds, "june")

summary(outer_results$roc_auc)

n_perms <- 5
perm_list <- rep(list(data), n_perms)

## FUTUREize this and other "map" functions, also abstract this out into a "run_permutations" function or something 
## (and make the other one "run_single_perm")
# perm_means <- future_map_dbl(perm_list, run_single_permutation, num_folds = 5, num_repeats = 5, desired_outcome = "20230824", 
#                              .options = furrr_options(seed = TRUE))

perm_means <- map_vec(perm_list, run_single_permutation, num_folds = 5, num_repeats = 5, desired_outcome = "20230824")

## Write perm_means to disk
perm_means %>% readr::write_rds(paste0(results_dir, "ni_lasso_date_perms_n", n_perms, "_", format(Sys.Date(), format="%Y%m%d"), ".rds"))

perm_means

hist(perm_means)
p_val <- sum (perm_means >= mean(outer_results$roc_auc)) / length(perm_means)



## Plot the Inner Fold Results
library(scales)

pooled_inner <- tuning_results %>% bind_rows



p <- ggplot(pooled_inner, aes(x = penalty, y = mean_roc)) + 
  scale_x_continuous(trans = 'log2') +
  xlab("LASSO Penalty") + ylab("Inner ROC")

for (i in 1:length(tuning_results)) {
  p <- p  +
    geom_line(data = tuning_results[[i]], alpha = .2) +
    geom_point(data = best_penalty(tuning_results[[i]]), pch = 16, alpha = 3/4)
}

p <- p + geom_smooth(data = pooled_inner, method = "loess", se = TRUE)
ggsave(paste0(figures_dir, "lasso_immune_dates.pdf"), height = 5, width = 5)


###------------------------------------------------------------