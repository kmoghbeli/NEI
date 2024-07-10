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

# Parallell Processing
if(parallel::detectCores(logical=FALSE) > 3) {
  library(doParallel)
  
  num_cores <- parallel::detectCores(logical=FALSE)
  cl <- makePSOCKcluster(num_cores - 2)
  registerDoParallel(cl)
}

source("OmicsToolbox/marker_genes.R")

##########################################

# Read in Seurat Objects
#june_tg <- SeuratDisk::LoadH5Seurat(paste0(data_dir, "june_tg_sct.h5Seurat"))
#aug_tg <- SeuratDisk::LoadH5Seurat(paste0(data_dir, "aug_tg_sct.h5Seurat"))
#combined_tg <- SeuratDisk::LoadH5Seurat(paste0(data_dir, "combined_tg_annot.h5Seurat"))
combined_tg <- SeuratDisk::LoadH5Seurat(paste0(data_dir, "combined_tg.h5Seurat"))

## Factorizations
combined_tg$condition <- factor(combined_tg$condition, levels = c("control", "scratch", "kos", "re"))
combined_tg$id <- factor(combined_tg$id)
combined_tg$date <- factor(combined_tg$date)
combined_tg$location <- factor(combined_tg$location)
combined_tg$celltype <- factor(combined_tg$celltype)

##########################################

# Compare Dates in the Integrated File

source("OmicsToolbox/lasso_nested_cv.R")

# subset1 = sample(Cells(combined_tg)[which(combined_tg$date == "20230629" & 
#                                             combined_tg$seurat_clusters == 0 & 
#                                             combined_tg$condition == "re")], 
#                  size = 500)
# 
# subset2 = sample(Cells(combined_tg)[which(combined_tg$date == "20230824" & 
#                                             combined_tg$seurat_clusters == 0 & 
#                                             combined_tg$condition == "re")], 
#                  size = 500)

subset1 = sample(Cells(combined_tg)[which(combined_tg$date == "20230629")], size = 500)
subset2 = sample(Cells(combined_tg)[which(combined_tg$date == "20230824")], size = 500)

lasso_data_subset <- subset(combined_tg, 
                            cells = c(subset1, subset2), 
                            features = VariableFeatures(combined_tg)) 

data <- lasso_data_subset@assays$SCT$data %>% as.matrix() %>% t() %>% 
  as.data.frame() %>% 
  rownames_to_column("barcode") %>% 
  left_join(lasso_data_subset@meta.data %>% select(date) %>% rownames_to_column("barcode"), 
            join_by("barcode" == "barcode")) %>% 
  rename(y = date) %>% 
  mutate(y = factor(y)) %>% 
  column_to_rownames("barcode")


### Nested CV
###------------------------------------------------------------
data_folds <- nested_cv(data, 
                        outside = vfold_cv(v = 5, repeats = 5, strata = y), 
                        inside = vfold_cv(v = 5, strata = y))

nestedcv_folds <- data_folds
desired_outcome <- "20230824"


# Step 1) Execute all the inner resampling loops:
# The object tuning_results is a list of data frames for each of the 50 outer resamples.
tuning_results <- future_map(nestedcv_folds$inner_resamples, 
                             summarize_tune_results, 
                             desired_outcome = desired_outcome, 
                             .options = furrr_options(seed = TRUE))

# Step 2) Collect the best penalties identified from the inner folds/resamples
best_penalty <- function(results) results[which.max(results$mean_roc),]

penalty_vals <- tuning_results %>% map_df(best_penalty) %>% select(penalty)

folds_with_penalty <- bind_cols(nestedcv_folds, penalty_vals)

# Step 3) Now that we have these estimates, we can compute the outer resampling results 
# for each of the outer splits using the corresponding tuning parameter value:

outer_results <- folds_with_penalty %>% 
  mutate(roc_auc = map2_dbl(splits, penalty, lasso_roc, desired_outcome = desired_outcome))

summary(outer_results$roc_auc)

################################################

### Lasso
###------------------------------------------------------------
X <- model.matrix(y ~ ., data = data)[,-1]
Y <- data[, "y"]

lasso.model <- glmnet::glmnet(x=X, y=Y, alpha = 1, lambda = 0.02, family = "binomial")

lasso_features <- lasso.model$beta %>% as_tibble(rownames = "feature") %>% 
  filter(s0 != 0) %>% 
  rename(coeff = s0) %>% 
  mutate(abs_coeff = abs(coeff)) %>% 
  arrange(desc(abs_coeff))

n <- 10
VlnPlot(combined_tg, assay = "SCT", 
        features = lasso_features$feature[1:n], 
        group.by = "date",
        scCustomize::DiscretePalette_scCustomize(n, palette = "glasbey", shuffle = FALSE),
        stack = TRUE, flip = TRUE) + NoLegend() + 
  ggtitle(paste0("Top ", n, " Lasso Features"))


#date_degs <- FindMarkers(combined_tg, assay = "SCT", group.by = "date", ident.1 = "20230629", ident.2 = "20230824", recorrect_umi = FALSE)

###------------------------------------------------------------


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

p + geom_smooth(data = pooled_inner, method = "loess", se = TRUE)
#ggsave(paste0(figures_dir, "lasso_immune_dates.pdf"), height = 5, width = 5)

# Step 4) Run Permutations (if desired)
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



##########################################

marker_genes_tg <- 
  c("Ptprc", 
    "Pax5", "Cd19", "Ighg1", "Ighm", "Ighd", "Cd27", # B cells - Naive (mouse) B cell markers (IgD+, CD27-)
    "Cd3d", "Cd3e", "Cd3g", # T
    "Klrb1c", "Prf1", "Klrk1", "Gzma", "Gzmb",  # NK 
    "Itga2", "Ncam1",  #NK-T
    "Itgam", "Adgre1", "Itgax",
    "Itgam", # Cd11b - Myeloid marker (mouse)
    "Adgre1", "Cd83", # F4/80, 
    "Itgax", # DCs
    "Cd14", "Cd68", "Fcgr1",  # Macs - note that Cd16 never comes up 
    "Ly6c1", 
    #"Cd74", # MHC-II mouse marker (used by Renthal 2022 to identify immune cells in TG)
    "Ptgs2", "Irf5", "Nos2",  # Mouse M1 Mac Markers 
    # "Stat1", "Retnla",  # Mouse M1 Mac Markers (less helpful)
    #"Il12a", "Il23a", "Cd163",  # M1 vs M2 (M1: IL-12 and IL23 high with CD163 neg and M2 the opposite)
    "Cd163",  # M2
    #"Arg1", # M2a
    "Socs3", "Cd86", # M2b
    "Ccr2", "Slamf6",   #M2c
    # "Tlr1", "Tlr8", "Scarb1", #M2c (less helpful)
    "Vegfa",    # M2d, 
    "Cx3cr1"  # Tissue-res Mac
  )

DimPlot(combined_tg, group.by = "seurat_clusters", 
        reduction = "umap", 
        #label = TRUE, label.box = TRUE, label.size = 3, repel = TRUE, 
        #cols = DiscretePalette(length(unique(combined_tg$seurat_clusters)), palette = "glasbey")
) + 
  NoLegend() + ggtitle("Clusters")
ggsave(paste0(figures_dir, "combined_tg_umap_seurat_clusters.png"), height = 7, width = 7)

VlnPlot(combined_tg, assay = "SCT",
        features = marker_genes_tg, 
        group.by = "seurat_clusters",
        stack = TRUE, flip = TRUE) + NoLegend() +
  ggtitle("Combined TG")

FeaturePlot(combined_tg, 
            features = c("Pax5", "Cd19", "Cd20", "Ighg1", "Ighm", "Ighd"), 
            pt.size = 1.4)

##########################################

# combined_tg_avgs <- AverageExpression(combined_tg, assays = "SCT", layer = "scale.data", group.by = "seurat_clusters",
#                                       features = c("Cd3d", "Klrb1c", "Pax5", "Cd68", "Nos2", "Irf5", "Cd163", "Socs3", "Ccr2", "Cx3cr1"))[[1]] %>%
#   as.matrix() %>% t() %>%
#   as_tibble(rownames = "cluster") %>%
#   mutate(designation = case_when(Pax5 > 0.3 ~ "B",
#                                  Cd3d > 1 ~ "T",
#                                  Klrb1c > 1 ~ "NK",
#                                  Cd68 > 0 & Nos2 > 0  & Cx3cr1 > 0 ~ "TR M1-like Mac",
#                                  Cd68 > 0 & Nos2 > 0 ~ "M1-like Mac",
#                                  Cd163 > 0 & Cx3cr1 > 0  ~ "TR M2a-like Mac",
#                                  Cd163 > 0 ~ "M2a-like Mac",
#                                  Cd68 > 0 & Socs3 > 0 & Cx3cr1 > 0   ~ "TR M2b-like Mac",
#                                  Cd68 > 0 & Socs3 > 0 ~ "M2b-like Mac",
#                                  Irf5 > 0 ~ "M1-like Mac",
#                                  .default = "other"
#                                  #.default = cluster
#                                  ),
#          .after = cluster)

combined_tg_avgs <- AverageExpression(combined_tg, assays = "SCT", layer = "scale.data", group.by = "seurat_clusters",
                                      features = c("Cd3d", "Klrb1c", "Pax5", "Cd68"))[[1]] %>%
  as.matrix() %>% t() %>%
  as_tibble(rownames = "cluster") %>%
  mutate(designation = case_when(Pax5 > 0.3 ~ "B",
                                 Cd3d > 1 ~ "T",
                                 Klrb1c > 1 ~ "NK",
                                 #Cd68 > 0 ~ "Mye",
                                 .default = "Mye"
                                 #.default = cluster
  ),
  .after = cluster)

print(combined_tg_avgs, n = Inf)

combined_tg$cell_L1 <- factor("immune")

combined_tg$cell_L2 <- factor(combined_tg_avgs$designation[as.numeric(combined_tg$seurat_clusters)], 
                              levels = c("Mye", "NK", "B", "T"))

# Cell_L4 (e.g., "Mac_0", "Mac_1")
combined_tg_cell_L4_mappings <- 
  combined_tg@meta.data %>% 
  select(cell_L2, seurat_clusters) %>% 
  arrange(cell_L2, seurat_clusters) %>% 
  group_by(cell_L2) %>% 
  distinct() %>% 
  mutate(idx = row_number() - 1, 
         cell_L4 = paste0(cell_L2, "_", idx)) %>% 
  select(-idx) %>% 
  right_join(combined_tg@meta.data %>% select(cell_L2, seurat_clusters) %>% rownames_to_column("barcode"), 
             join_by("cell_L2" == "cell_L2", "seurat_clusters" == "seurat_clusters")) %>% 
  column_to_rownames("barcode") %>% 
  select(cell_L4)

combined_tg <- AddMetaData(combined_tg, combined_tg_cell_L4_mappings)
combined_tg$cell_L4 <- factor(combined_tg$cell_L4, 
                              levels =  str_sort(unique(combined_tg$cell_L4), numeric = TRUE))

## Take a quick view at the violin plot now to see if things line up as they should
VlnPlot(combined_tg, assay = "SCT",
        features = marker_genes_tg, 
        group.by = "cell_L4",
        stack = TRUE, flip = TRUE) + NoLegend() +
  ggtitle("Combined TG")


## Write the new annotated version and delete our unannotated version
combined_tg %>% SeuratDisk::SaveH5Seurat(paste0(data_dir, "combined_tg_annot.h5Seurat"), overwrite = TRUE)

file.remove(paste0(data_dir, "combined_tg.h5Seurat"))  # Delete the old file

##########################################

# ggsave(paste0(figures_dir, "combined_tg_markers.png"), height = 7, width = 7)

VlnPlot(combined_tg, assay = "SCT",
        features = marker_genes_tg,
        group.by = "cell_L2",
        scCustomize::DiscretePalette_scCustomize(length(marker_genes_tg), palette = "glasbey", shuffle = FALSE),
        stack = TRUE, flip = TRUE) + NoLegend() +
  ggtitle("Combined TG")
# 
# ggsave(paste0(figures_dir, "combined_tg_markers_with_annot.png"), height = 7, width = 7)

VlnPlot(combined_tg, assay = "SCT",
        features = marker_genes_tg,
        group.by = "cell_L3",
        scCustomize::DiscretePalette_scCustomize(length(marker_genes_tg), palette = "glasbey", shuffle = FALSE),
        stack = TRUE, flip = TRUE) + NoLegend() +
  ggtitle("Combined TG")


# Idents(combined_tg) <- combined_tg$cell_L2
# not_other <- levels(combined_tg$cell_L2)["other" != levels(combined_tg$cell_L2)]
# combined_tg_no_other <- subset(combined_tg, idents = not_other)
# 
# DimPlot(combined_tg_no_other, group.by = "cell_L2", reduction = "umap", label = TRUE, label.box = TRUE, label.size = 3, repel = TRUE) + 
#   ggtitle("Trigeminal Ganglion") + NoLegend()
# 
# ggsave(paste0(figures_dir, "combined_tg_umap.png"), height = 7, width = 7)


# Sanity check for all of these (should all equal to 100)
# combined_tg@meta.data %>% 
#   group_by(condition, cell_L2) %>% 
#   summarise(n = n(), .groups = "drop_last") %>%
#   mutate(prop = n / sum(n) * 100) %>% 
#   summarise(summm = sum(prop))

# Cell L2 and L3 by condition (just for all immune cells)
combined_tg@meta.data %>% 
  group_by(condition, cell_L2) %>% 
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(prop = n / sum(n) * 100) %>% 
  ggplot(aes(cell_L2, prop, fill = condition)) + 
  scale_fill_viridis_d() + 
  geom_col(position=position_dodge()) + 
  ggprism::theme_prism() + 
  labs(x = "", y = "Proportion of all immune cells", title = "Trigeminal Ganglion")

ggsave(paste0(figures_dir, "combined_tg_props_immune_L2.png"), height = 7, width = 14)

# Control Only
combined_tg@meta.data %>% 
  filter("control" == condition) %>% 
  group_by(condition, cell_L2) %>% 
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(prop = n / sum(n) * 100) %>% 
  ggplot(aes(cell_L2, prop)) + 
  scale_fill_viridis_d() + 
  geom_col(position=position_dodge()) + 
  ggprism::theme_prism() + 
  labs(x = "", y = "Proportion of TG immune cells") + 
  theme(aspect.ratio = 1.7)

ggsave(paste0(figures_dir, "combined_tg_props_immune_L2_control_only.png"), height = 7, width = 14)

# Control only by date
combined_tg@meta.data %>% 
  filter("control" == condition) %>% 
  group_by(date, condition, cell_L2) %>% 
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(prop = n / sum(n) * 100) %>% 
  ggplot(aes(cell_L2, prop, fill = date)) + 
  geom_col(position=position_dodge()) + 
  ggprism::theme_prism() + 
  ggprism::scale_fill_prism(palette = "colorblind_safe") + 
  labs(x = "", y = "Proportion of TG immune cells") + 
  theme(aspect.ratio = 1.7)

# Fisher's Exact Test - Cell type proportions in Control by date
control_metadata <- combined_tg@meta.data %>% filter("control" == condition)

fisher.test(t(table(control_metadata$date, control_metadata$cell_L2)), workspace = 200000000)



combined_tg@meta.data %>% 
  # filter("immune" == cell_L1, 
  #        condition %in% c("kos", "re")) %>% 
  group_by(date, condition, cell_L2) %>% 
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(prop = n / sum(n) * 100) %>%
  ggplot(aes(cell_L2, prop, fill = condition)) + 
  scale_fill_viridis_d() + 
  geom_col(width = 0.55, position=position_dodge()) + 
  facet_grid(~ date) + 
  ggprism::theme_prism() + 
  #theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  labs(x = "", y = "Proportion of all immune cells", title = "Trigeminal Ganglion")

ggsave(paste0(figures_dir, "combined_tg_props_immune_L2_date.png"), height = 7, width = 14)


combined_tg@meta.data %>% 
  group_by(condition, cell_L3) %>% 
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(prop = n / sum(n) * 100) %>% 
  ggplot(aes(cell_L3, prop, fill = condition)) + 
  scale_fill_viridis_d() + 
  geom_col(position=position_dodge()) + 
  ggprism::theme_prism() + 
  labs(x = "", y = "Proportion of all immune cells", title = "Trigeminal Ganglion")

ggsave(paste0(figures_dir, "combined_tg_props_immune_L3.png"), height = 7, width = 14)

combined_tg@meta.data %>% 
  # filter("immune" == cell_L1, 
  #        condition %in% c("kos", "re")) %>% 
  group_by(date, condition, cell_L3) %>% 
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(prop = n / sum(n) * 100) %>%
  ggplot(aes(cell_L3, prop, fill = condition)) + 
  scale_fill_viridis_d() + 
  geom_col(width = 0.55, position=position_dodge()) + 
  facet_grid(~ date) + 
  ggprism::theme_prism() + 
  #theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  labs(x = "", y = "Proportion of all immune cells", title = "Trigeminal Ganglion")

ggsave(paste0(figures_dir, "combined_tg_props_immune_L3_date.png"), height = 7, width = 20)


combined_tg@meta.data %>% 
  filter("immune" == cell_L1, 
          condition %in% c("kos", "re")) %>% 
  group_by(date, condition, cell_L3) %>% 
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(prop = n / sum(n) * 100) %>%
  ggplot(aes(cell_L3, prop, fill = condition)) + 
  scale_fill_viridis_d(end = 0.75) + 
  geom_col(width = 0.55, position=position_dodge()) + 
  facet_grid(~ date) + 
  ggprism::theme_prism() + 
  #theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  labs(x = "", y = "Proportion of all immune cells", title = "Trigeminal Ganglion")

ggsave(paste0(figures_dir, "combined_tg_props_immune_L3_date_re_kos.png"), height = 7, width = 21)


table(combined_tg$seurat_clusters, combined_tg$condition)

##########################################

Idents(combined_tg) <- combined_tg$condition
combined_tg_no_scratch <- subset(combined_tg, idents = c("control", "kos", "re"))


brian_genes <- c("Ccr2", 
                 "Ifngr1", "Ifngr2", "Ifnar1", "Ifnar2", 
                 "Tnfaip1", "Tnfaip3", "Tnfrsf21", "Tnfrsf11a", "Tnfrsf1a", 
                 "Il4ra", "Il6", "Il6st", "Il10rb", "Il13ra1"
)

VlnPlot(combined_tg_no_scratch, assay = "SCT",
        features = brian_genes,
        group.by = "cell_L3", split.by = "condition", 
        cols = viridis::viridis(4),
        stack = TRUE, flip = TRUE) + 
  ggtitle("Combined TG") #+ ggprism::theme_prism()

ggsave(paste0(figures_dir, "combined_tg_brian_genes.png"), height = 7, width = 14)


nt_receptor_genes <- c("Nmur1", "Adrb2", 
                       #"Calcr", 
                       "Vipr2", "Chrm1", 
                       #"Chrna1", "Chrna7", 
                       "Chrnb1", "Fpr1", 
                       #"Mrgpra1", 
                       "Mrgprb1", "Grm5", "Grm7", "Grin1", "Grin2a", "Tacr1", 
                       "Bmpr2", "Csf3r")

VlnPlot(combined_tg_no_scratch, assay = "SCT",
        features = nt_receptor_genes,
        group.by = "cell_L2", split.by = "condition", 
        cols = viridis::viridis(4),
        stack = TRUE, flip = TRUE) + 
  ggtitle("Combined TG") #+ ggprism::theme_prism()

ggsave(paste0(figures_dir, "combined_tg_ntr_genes_L2.png"), height = 5, width = 5)

VlnPlot(combined_tg_no_scratch, assay = "SCT",
        features = nt_receptor_genes,
        group.by = "cell_L3", split.by = "condition", 
        cols = viridis::viridis(4),
        stack = TRUE, flip = TRUE) + 
  ggtitle("Combined TG") #+ ggprism::theme_prism()

ggsave(paste0(figures_dir, "combined_tg_ntr_genes_L3.png"), height = 7, width = 10)

VlnPlot(combined_tg_no_scratch, assay = "SCT",
        features = c("Sarm1"),
        group.by = "cell_L2", split.by = "condition", 
        cols = viridis::viridis(4),
        flip = TRUE) + 
  ggtitle("Combined TG") #+ ggprism::theme_prism()

ggsave(paste0(figures_dir, "combined_tg_sarm1_L2.png"), height = 5, width = 5)

##########################################



