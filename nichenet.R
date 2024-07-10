library(tidyverse)
library(ggprism)
library(Seurat)
library(nichenetr)
library(patchwork)
library(ComplexHeatmap)
options(stringsAsFactors = FALSE)
conflicted::conflict_prefer_all("dplyr", quiet = TRUE)

kaveh_colors1 <- c("#a6611a", "#dfc27d", "#f5f5f5", "#80cdc1", "#018571")
kaveh_colors2 <- c("#40004b", "#762a83", "#9970ab", "#c2a5cf", "#e7d4e8", "#f7f7f7", "#d9f0d3", "#a6dba0", "#5aae61", "#1b7837", "#00441b")
kaveh_colors3 <- c("#ff0000", "#f5f5f5","#0000ff")


DATASETS_DIR <- "~/Library/CloudStorage/Dropbox/Projects/_DATASETS/NicheNetV2/"

config_dir <- "./config/"
data_dir <- "./data/"
model_dir <- "./models/"
results_dir <- "./results/"
figures_dir <- "./figures/"

## CUSTOM GGPLOT THEME
theme_set(ggprism::theme_prism(palette = "colorblind_safe", 
                               base_size = 16, 
                               base_line_size = 1.5))

# Parallell Processing
if(parallel::detectCores(logical=FALSE) > 3) {
  library(doParallel)
  
  num_cores <- parallel::detectCores(logical=FALSE)
  cl <- makePSOCKcluster(num_cores - 2)
  registerDoParallel(cl)
}

#######################################################################################################################################

dataset <- "sept_cornea"
#dataset <- "combined_tg"
grouping <- "cell_L2"
#grouping <- "cell_L3"
#grouping <- "cell_L4"

plot_prefix <- paste0("nn_", dataset, "_", grouping, "_")

## Set Sender / Receiver groups of Interest
if (grepl("tg", dataset, ignore.case = TRUE)) {
  
  # Tg: Mye -> Mye
  
  receiver = "Mye"
  sender_celltypes <- c("Mye")
  
} else if(grepl("cornea", dataset, ignore.case = TRUE)) {
  
  # Cornea: NK -> Mye, B -> Mye, B -> Epi, T->T, T->Epi
  
  receiver = "Mye"
  sender_celltypes <- c("NK", "Mye")
}

## https://github.com/saeyslab/nichenetr/blob/master/vignettes/seurat_steps.md
## 1) Load the full annotated, normalized dataset

seurat_data <- SeuratDisk::LoadH5Seurat(paste0(data_dir, dataset, "_annot.h5Seurat"))  
Idents(seurat_data) <- seurat_data[[grouping]] %>% rownames_to_column() %>% deframe()

## Read in NicheNet’s networks
# The ligand-target prior model is a matrix describing the potential that a ligand may regulate a target gene, and it is used to run the ligand activity analysis. 
# The ligand-receptor network contains information on potential ligand-receptor bindings, and it is used to identify potential ligands. 
# The weighted ligand-receptor network contains weights representing the potential that a ligand will bind to a receptor, and it is used for visualization.

organism <- "mouse"

if(organism == "human") {
  lr_network <- readRDS(file = paste0(DATASETS_DIR, "lr_network_human_21122021.rds"))
  ligand_target_matrix <- readRDS(file = paste0(DATASETS_DIR, "ligand_target_matrix_nsga2r_final.rds"))
  weighted_networks <- readRDS(file = paste0(DATASETS_DIR, "weighted_networks_nsga2r_final.rds"))
  
  ligand_tf_matrix <- readRDS(file = paste0(DATASETS_DIR, "ligand_tf_matrix_nsga2r_final.rds"))
  sig_network <- readRDS(file = paste0(DATASETS_DIR, "signaling_network_human_21122021.rds"))
  gr_network <- readRDS(file = paste0(DATASETS_DIR, "gr_network_human_21122021.rds"))
  
} else if(organism == "mouse") {
  lr_network <- readRDS(file = paste0(DATASETS_DIR, "lr_network_mouse_21122021.rds"))
  ligand_target_matrix <- readRDS(file = paste0(DATASETS_DIR, "ligand_target_matrix_nsga2r_final_mouse.rds"))
  weighted_networks <- readRDS(file = paste0(DATASETS_DIR, "weighted_networks_nsga2r_final_mouse.rds"))
  
  ligand_tf_matrix <- readRDS(file = paste0(DATASETS_DIR, "ligand_tf_matrix_nsga2r_final_mouse.rds"))
  sig_network <- readRDS(file = paste0(DATASETS_DIR, "signaling_network_mouse_21122021.rds"))
  gr_network <- readRDS(file = paste0(DATASETS_DIR, "gr_network_mouse_21122021.rds"))
}


## Perform the NicheNet analysis
expressed_genes_receiver <- get_expressed_genes(receiver, seurat_data, pct = 0.05)

lr_network <- lr_network %>% distinct(from, to)

## 1. Define a set of potential ligands for both the sender-agnostic and sender-focused approach

# Get a list of all receptors available in the ligand-receptor network, 
# and define expressed receptors as genes that are in the ligand-receptor 
# network and expressed in the receiver. Then, define the potential ligands 
# as all ligands whose cognate receptors are expressed.

all_receptors <- unique(lr_network$to)
expressed_receptors <- intersect(all_receptors, expressed_genes_receiver)
potential_ligands <- lr_network %>% filter(to %in% expressed_receptors) %>% pull(from) %>% unique()

# Sender-focused approach
# Use lapply to get the expressed genes of every sender cell type separately here
list_expressed_genes_sender <- sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seurat_data, 0.05)
expressed_genes_sender <- list_expressed_genes_sender %>% unlist() %>% unique()

potential_ligands_focused <- intersect(potential_ligands, expressed_genes_sender) 

# Also check 
length(expressed_genes_sender)
length(potential_ligands)
length(potential_ligands_focused)


## 2. Define the gene set of interest
metadata_var <- "condition"
condition_oi <-  "re"
condition_reference <- "kos"

seurat_obj_receiver <- subset(seurat_data, idents = receiver)

DE_table_receiver <-  FindMarkers(object = seurat_obj_receiver,
                                  ident.1 = condition_oi, ident.2 = condition_reference,
                                  group.by = metadata_var, 
                                  recorrect_umi = FALSE,
                                  min.pct = 0.05) %>% rownames_to_column("gene")

geneset_oi <- DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
geneset_oi <- geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]


## 3. Define the background genes
# All expressed genes in the receiver cell population (that are also in the 
# ligand-target matrix) is defined as the ‘background set’ for our ligand 
# prioritization procedure in the next step. It’s also important to check that 
# the number of background genes is a ‘reasonable’ number, generally between 5000-10000, 
# and sufficiently larger than our gene set of interest
background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

length(background_expressed_genes)
length(geneset_oi)

## 4. Perform NicheNet ligand activity analysis
# This is the main step of NicheNet where the potential ligands are ranked 
# based on the presence of their target genes in the gene set of interest 
# (compared to the background set of genes). 
# Ligands are ranked based on the area under the precision-recall curve (AUPR) between a ligand’s target predictions and the observed transcriptional response. Although other metrics like the AUROC and pearson correlation coefficient are also computed, we demonstrated in our validation study that the AUPR was the most informative measure to define ligand activity (this was the Pearson correlation for v1).

ligand_activities <- predict_ligand_activities(geneset = geneset_oi,
                                               background_expressed_genes = background_expressed_genes,
                                               ligand_target_matrix = ligand_target_matrix,
                                               potential_ligands = potential_ligands)

ligand_activities <- ligand_activities %>% arrange(-aupr_corrected) %>% mutate(rank = rank(desc(aupr_corrected)))
ligand_activities

num_top_ligands <- 50
#num_top_ligands <- 100

# We will use the top 30 ligands to predict active target genes and construct 
# an active ligand-receptor network. However, the choice of looking only at 
# the 30 top-ranked ligands for further biological interpretation is based 
# on biological intuition and is quite arbitrary. Therefore, users can decide 
# to continue the analysis with a different number of ligands. We recommend 
# to check the selected cutoff by looking at the distribution of the ligand 
# activity values. Here, we show the ligand activity histogram (the score for 
# the 30th ligand is indicated via the dashed line).
ggplot(ligand_activities, aes(x=aupr_corrected)) + 
  geom_histogram(color="black", fill="darkorange")  + 
  geom_vline(aes(xintercept=min(ligand_activities %>% top_n(num_top_ligands, aupr_corrected) %>% pull(aupr_corrected))),
             color="red", linetype="dashed", linewidth=1) + 
  labs(x="ligand activity (PCC)", y = "# ligands") +
  theme_classic()

best_upstream_ligands <- ligand_activities %>% top_n(num_top_ligands, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand)

vis_ligand_aupr <- ligand_activities %>% filter(test_ligand %in% best_upstream_ligands) %>%
  column_to_rownames("test_ligand") %>% select(aupr_corrected) %>% arrange(aupr_corrected) %>% as.matrix(ncol = 1)

(make_heatmap_ggplot(vis_ligand_aupr,
                     "Prioritized ligands", "Ligand activity", 
                     legend_title = "AUPR", color = "darkorange") + 
    theme(axis.text.x.top = element_blank())) 

# 5. Infer target genes and receptors of top-ranked ligands
# Active target gene inference
# Active target genes are defined as genes in the gene set of interest that have the highest regulatory potential for each top-ranked ligand. 
# These top targets of each ligand are based on the prior model. 
# The function get_weighted_ligand_target_links will return genes that are in the gene set of interest and are the top n targets of a ligand (default: n = 200, but there are too many target genes here so we only considered the top 100).
active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,
         geneset = geneset_oi,
         ligand_target_matrix = ligand_target_matrix,
         n = 100) %>%
  bind_rows() %>% drop_na()

# For visualization purposes, the ligand-target prior model was adapted by 
# setting a regulatory potential score to 0 if their score was below a 
# predefined cutoff (default: 0.25, or the 25th percentile) across all 
# scores between the top-ranked ligands and their top n targets. We recommend 
# users to test several cutoff values for the best visualization, as lowering 
# or increasing the cutoff will result in a denser or sparser heatmap, respectively.

active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.85) 

order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))

vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])

make_heatmap_ggplot(vis_ligand_target, "Prioritized ligands", "Predicted target genes",
                    color = "purple", legend_title = "Regulatory\npotential\n") +
  scale_fill_gradient2(low = "whitesmoke",  high = "purple")

# Receptors of top-ranked ligands
# Similar to above, we identify which receptors have the highest interaction potential with the top-ranked ligands
ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
  best_upstream_ligands, expressed_receptors,
  lr_network, weighted_networks$lr_sig) 

vis_ligand_receptor_network <- prepare_ligand_receptor_visualization(
  ligand_receptor_links_df,
  best_upstream_ligands,
  order_hclust = "both") 

(make_heatmap_ggplot(t(vis_ligand_receptor_network), 
                     y_name = "Ligands", x_name = "Receptors",  
                     color = "mediumvioletred", legend_title = "Prior interaction potential"))

#ggsave(paste0(figures_dir, plot_prefix, "_predicted_receptors.png") , plot = last_plot(), width = 14, height = 7)


## 6. Sender-focused approach
# To perform the sender-focused approach, simply subset the ligand activities 
# to only contain expressed ligands from all populations (calculated in Step 1). 
# We can then perform target gene and receptor inference as above.
ligand_activities_all <- ligand_activities 
best_upstream_ligands_all <- best_upstream_ligands

ligand_activities <- ligand_activities %>% filter(test_ligand %in% potential_ligands_focused)
best_upstream_ligands <- ligand_activities %>% top_n(30, aupr_corrected) %>% arrange(-aupr_corrected) %>%
  pull(test_ligand) %>% unique()


ligand_aupr_matrix <- ligand_activities %>% filter(test_ligand %in% best_upstream_ligands) %>%
  column_to_rownames("test_ligand") %>% select(aupr_corrected) %>% arrange(aupr_corrected)

vis_ligand_aupr <- as.matrix(ligand_aupr_matrix, ncol = 1) 

p_ligand_aupr <- make_heatmap_ggplot(vis_ligand_aupr,
                                     "Prioritized ligands", "Ligand activity", 
                                     legend_title = "AUPR", color = "darkorange") + 
  theme(axis.text.x.top = element_blank())

p_ligand_aupr 
#ggsave(paste0(figures_dir, plot_prefix, "_ligand_activ.png") , plot = last_plot(), width = 5, height = 7)

# Target gene plot
active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,
         geneset = geneset_oi,
         ligand_target_matrix = ligand_target_matrix,
         n = 100) %>%
  bind_rows() %>% drop_na()

active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.90) 

order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))

vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])

p_ligand_target <- make_heatmap_ggplot(vis_ligand_target, "Prioritized ligands", "Predicted target genes",
                                       color = "purple", legend_title = "Regulatory\npotential\n") +
  scale_fill_gradient2(low = "whitesmoke",  high = "purple")

#ggsave(paste0(figures_dir, plot_prefix, "_predicted_target_genes.png") , plot = p_ligand_target, width = 14, height = 7)


# Receptor plot
ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
  best_upstream_ligands, expressed_receptors,
  lr_network, weighted_networks$lr_sig) 

vis_ligand_receptor_network <- prepare_ligand_receptor_visualization(
  ligand_receptor_links_df,
  best_upstream_ligands,
  order_hclust = "both") 

p_ligand_receptor <- make_heatmap_ggplot(t(vis_ligand_receptor_network), 
                                         y_name = "Ligands", x_name = "Receptors",  
                                         color = "mediumvioletred", legend_title = "Prior interaction potential")

#ggsave(paste0(figures_dir, plot_prefix, "_predicted_receptors.png") , plot = p_ligand_receptor, width = 14, height = 7)

# Visualizing expression and log-fold change in sender cells
# For the sender-focused approach, we can also investigate further on which 
# sender cell populations are potentially the true sender of these ligands. 
# First, we can simply check which sender cell population expresses which of 
# these top-ranked ligands.
# Dotplot of sender-focused approach
p_dotplot <- DotPlot(subset(seurat_data, cell_L2 %in% sender_celltypes),
                     features = rev(best_upstream_ligands), cols = "RdYlBu") + 
  coord_flip() + 
  xlab("Top-ranked Ligands") + 
  ylab("Cell Type") + 
  scale_y_discrete(position = "right")

p_dotplot
#ggsave(paste0(figures_dir, plot_prefix, "_potential_senders.png") , plot = last_plot(), width = 7, height = 7)


# Next, we can also check upregulation of ligands in sender cells by computing the log-fold change between the two conditions. This ligand differential expression is not used for prioritization and ranking of the ligands (the ranking is only determined based on enrichment of target genes among DE genes in the receiver, CD8T cells), but it can add a useful extra layer of information next to the ligand activities. This is of course only possible in some cases, such as case-control studies.
celltype_order <- levels(Idents(seurat_data)) 

# Use this if cell type labels are the identities of your Seurat object
# if not: indicate the celltype_col properly
DE_table_top_ligands <- lapply(
  celltype_order[celltype_order %in% sender_celltypes],
  get_lfc_celltype, 
  seurat_obj = seurat_data,
  condition_colname = metadata_var,
  condition_oi = condition_oi,
  condition_reference = condition_reference,
  celltype_col = grouping,
  min.pct = 0, logfc.threshold = 0, 
  recorrect_umi = FALSE,
  features = best_upstream_ligands 
) 

DE_table_top_ligands <- DE_table_top_ligands %>%  reduce(., full_join) %>% 
  column_to_rownames("gene") 

vis_ligand_lfc <- as.matrix(DE_table_top_ligands[rev(best_upstream_ligands), ]) 

p_lfc <- make_threecolor_heatmap_ggplot(vis_ligand_lfc,
                                        "Prioritized ligands", "LFC in Sender",
                                        low_color = "midnightblue", mid_color = "white",
                                        mid = median(vis_ligand_lfc), high_color = "red",
                                        legend_title = "LFC")

#ggsave(paste0(figures_dir, plot_prefix, "_sender_lfc.png") , plot = last_plot(), width = 7, height = 7)


# Finally, you can also compare rankings between the sender-agnostic and 
# sender-focused approach. Here, the red sections of the left bar plot indicates 
# which ligands in the sender-agnostic approach are filtered out in the 
# sender-focused approach because they are not expressed.

#### DOES NOT WORK
# (make_line_plot(ligand_activities = ligand_activities_all,
#                 potential_ligands = potential_ligands_focused) +
#     theme(plot.title = element_text(size=11, hjust=0.1, margin=margin(0, 0, -5, 0))))


## 7. Summary visualizations of the NicheNet analysis
# Finally, we can make a combined plot containing heatmap of ligand activities, ligand expression, ligand log-fold change and the target genes of the top-ranked ligands. As mentioned earlier, sometimes ligands do not appear in the ligand-target heatmap because they don’t have target genes with high enough regulatory potential scores. In this case, CCl22 is present in other plots (ranked 25th) but is missing in the rightmost plot. If users wish for these plots to be consistent, they may use the variable order_ligands defined when creating the ligand-target heatmap to subset other plots instead of best_upstream_ligands.
figures_without_legend <- cowplot::plot_grid(
  p_ligand_aupr + theme(legend.position = "none"),
  p_dotplot + theme(legend.position = "none",
                    axis.ticks = element_blank(),
                    axis.title.y = element_blank(),
                    axis.title.x = element_text(size = 12),
                    axis.text.y = element_text(size = 9),
                    axis.text.x = element_text(size = 9,  angle = 90, hjust = 0)) +
    ylab("Expression in Sender"),
  p_lfc + theme(legend.position = "none",
                axis.title.y = element_blank()),
  p_ligand_target + theme(legend.position = "none",
                          axis.title.y = element_blank()),
  align = "hv",
  nrow = 1,
  rel_widths = c(ncol(vis_ligand_aupr)+6, ncol(vis_ligand_lfc)+7, ncol(vis_ligand_lfc)+8, ncol(vis_ligand_target)))

legends <- cowplot::plot_grid(
  ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_aupr)),
  ggpubr::as_ggplot(ggpubr::get_legend(p_dotplot)),
  ggpubr::as_ggplot(ggpubr::get_legend(p_lfc)),
  ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_target)),
  nrow = 1,
  align = "h", rel_widths = c(1.5, 1, 1, 1))

combined_plot <-  cowplot::plot_grid(figures_without_legend, legends, rel_heights = c(10,5), nrow = 2, align = "hv")
combined_plot

ggsave(paste0(figures_dir, plot_prefix, "_summary.png") , plot = combined_plot, width = 25, height = 14, bg = "white")

## 8. Inferring ligand-to-target signaling paths
# https://github.com/saeyslab/nichenetr/blob/master/vignettes/ligand_target_signaling_path.md
# To determine signaling paths between a ligand and target of interest, we look at which transcription factors are best regulating the target genes and are most closely downstream of the ligand (based on the weights of the edges in the integrated ligand-signaling and gene regulatory networks). Then, the shortest paths between these transcription factors and the ligand of interests are determined and genes forming part in this path are considered as important signaling mediators. Finally, we look in our collected data source networks for all interactions between the ligand, signaling mediators, transcription factors and target genes. This allows to both prioritize signaling mediators and check which of all collected data sources support the ligand-target predictions of interest.

# Infer signaling paths between ligand and targets
#ligands_oi <- "Ifng" # this can be a list of multiple ligands if required
ligands_oi <- best_upstream_ligands[1] # this can be a list of multiple ligands if required
#targets_oi <- c("Stat1", "Irf1")
targets_oi <- active_ligand_target_links_df %>% 
  filter(ligand %in% ligands_oi) %>% 
  group_by(ligand) %>% 
  arrange(desc(weight)) %>% 
  slice_head(n = 3) %>% 
  pull(target) %>% 
  sort()

active_signaling_network <- get_ligand_signaling_path(ligands_all = ligands_oi,
                                                      targets_all = targets_oi,
                                                      weighted_networks = weighted_networks,
                                                      ligand_tf_matrix = ligand_tf_matrix,
                                                      top_n_regulators = 4,
                                                      minmax_scaling = TRUE) 


graph_min_max <- diagrammer_format_signaling_graph(signaling_graph_list = active_signaling_network,
                                                   ligands_all = ligands_oi, targets_all = targets_oi,
                                                   sig_color = "indianred", gr_color = "steelblue")

# To render the graph in RStudio Viewer, uncomment following line of code
# DiagrammeR::render_graph(graph_min_max, layout = "tree")

# To export/draw the svg, you need to install DiagrammeRsvg
graph_svg <- DiagrammeRsvg::export_svg(DiagrammeR::render_graph(graph_min_max, layout = "tree", output = "graph"))
interaction_plot <- cowplot::ggdraw() + cowplot::draw_image(charToRaw(graph_svg))

ggsave(paste0(figures_dir, plot_prefix, "_ligand_target_sig_path.png") , plot = interaction_plot, width = 7, height = 7, bg = "white")

# We will now look which of the collected data sources support the interactions in this network.
data_source_network <- infer_supporting_datasources(signaling_graph_list = active_signaling_network,
                                                    lr_network = lr_network, sig_network = sig_network, gr_network = gr_network)
head(data_source_network) 

# Export to Cytoscape: export the following to e.g. Cytoscape for exploration of the networks
# output_path <- ""
# write_output <- FALSE # change to TRUE for writing output
# 
# # weighted networks ('import network' in Cytoscape)
# if(write_output){
#   bind_rows(active_signaling_network$sig %>% mutate(layer = "signaling"),
#             active_signaling_network$gr %>% mutate(layer = "regulatory")) %>%
#     write_tsv(paste0(output_path,"weighted_signaling_network.txt")) 
# }
# 
# # networks with information of supporting data sources ('import network' in Cytoscape)
# if(write_output){
#   data_source_network %>% write_tsv(paste0(output_path,"data_source_network.txt"))
# }
# 
# # Node annotation table ('import table' in Cytoscape)
# specific_annotation_tbl <- bind_rows(
#   tibble(gene = ligands_oi, annotation = "ligand"),
#   tibble(gene = targets_oi, annotation = "target"),
#   tibble(gene = c(data_source_network$from, data_source_network$to) %>% unique() %>% setdiff(c(targets_oi,ligands_oi)) %>% intersect(lr_network$to %>% unique()), annotation = "receptor"),
#   tibble(gene = c(data_source_network$from, data_source_network$to) %>% unique() %>% setdiff(c(targets_oi,ligands_oi)) %>% intersect(gr_network$from %>% unique()) %>% setdiff(c(data_source_network$from, data_source_network$to) %>% unique() %>% intersect(lr_network$to %>% unique())),annotation = "transcriptional regulator")
# )
# non_specific_annotation_tbl <- tibble(gene = c(data_source_network$from, data_source_network$to) %>% unique() %>% setdiff(specific_annotation_tbl$gene), annotation = "signaling mediator")
# 
# if(write_output){
#   bind_rows(specific_annotation_tbl, non_specific_annotation_tbl) %>%
#     write_tsv(paste0(output_path,"annotation_table.txt"))
# }


## 9. Assess how well top-ranked ligands can predict a gene set of interest
# For the top 30 ligands, we will now build a multi-ligand model that uses all 
# top-ranked ligands to predict whether a gene belongs to the gene set of 
# interest (differentially expressed genes in CD8 T cells after LCMV infection) 
# or not. This classification model will be trained via cross-validation and 
# returns a probability for every gene.

# Change rounds and folds here, to two rounds to reduce time: normally: do multiple rounds
k <- 3 # 3-fold
n <- 2 # 2 rounds

gene_predictions_top30_list <- lapply(1:n, assess_rf_class_probabilities,
                                      folds = k,
                                      geneset = geneset_oi,
                                      background_expressed_genes = background_expressed_genes,
                                      ligands_oi = best_upstream_ligands,
                                      ligand_target_matrix = ligand_target_matrix)

# get performance: auroc-aupr-pearson
target_prediction_performances_cv <- gene_predictions_top30_list %>% lapply(classification_evaluation_continuous_pred_wrapper) %>%
  bind_rows() %>% mutate(round=seq(1:nrow(.)))

target_prediction_performances_cv$auroc %>% mean()
target_prediction_performances_cv$aupr %>% mean()
target_prediction_performances_cv$pearson %>% mean()

# get performance: how many viral response genes and non-viral response-genes among top 5% predicted targets
target_prediction_performances_discrete_cv <- gene_predictions_top30_list %>%
  lapply(calculate_fraction_top_predicted,
         quantile_cutoff = 0.95) %>%
  bind_rows(.id = "round") 

# What is the fraction of viral response genes that belongs to the top 5% predicted targets?
  
target_prediction_performances_discrete_cv %>% filter(true_target) %>% .$fraction_positive_predicted %>% mean()

# What is the fraction of non-viral-response genes that belongs to the top 5% predicted targets?
target_prediction_performances_discrete_cv %>% filter(!true_target) %>% .$fraction_positive_predicted %>% mean()

# We see that the viral response genes are enriched in the top-predicted target genes. To test this, we will now apply a Fisher’s exact test for every cross-validation round and report the average p-value.
target_prediction_performances_discrete_fisher <- gene_predictions_top30_list %>%
  lapply(calculate_fraction_top_predicted_fisher,
         quantile_cutoff = 0.95)

target_prediction_performances_discrete_fisher %>% unlist() %>% mean()

# get top predicted genes
top_predicted_genes <- lapply(1:n, get_top_predicted_genes,
                              gene_predictions_top30_list) %>%
  reduce(full_join, by = c("gene","true_target"))

top_predicted_genes %>% filter(true_target)

