library(tidyverse)
library(ggprism)
library(Seurat)
library(SeuratDisk)
library(CellChat)
library(patchwork)
library(ComplexHeatmap)
options(stringsAsFactors = FALSE)
conflicted::conflict_prefer_all("dplyr", quiet = TRUE)

kaveh_colors1 <- c("#a6611a", "#dfc27d", "#f5f5f5", "#80cdc1", "#018571")
kaveh_colors2 <- c("#40004b", "#762a83", "#9970ab", "#c2a5cf", "#e7d4e8", "#f7f7f7", "#d9f0d3", "#a6dba0", "#5aae61", "#1b7837", "#00441b")
kaveh_colors3 <- c("#ff0000", "#f5f5f5","#0000ff")


config_dir <- "./config/"
data_dir <- "../data_objects/"
figures_dir <- "../figures/"

## CUSTOM GGPLOT THEME
GG_KM_THEME <- 
  list(
    ggprism::theme_prism(palette = "colorblind_safe", 
                         base_size = 16, 
                         base_line_size = 1.5),
    scale_fill_viridis_d(alpha = 0.7, begin = 0.14, end = 0.88)
  )

# Parallell Processing
if(parallel::detectCores(logical=FALSE) > 3) {
  library(doParallel)
  
  num_cores <- parallel::detectCores(logical=FALSE)
  cl <- makePSOCKcluster(num_cores - 2)
  registerDoParallel(cl)
}

#######################################################################################################################################

dataset <- "tg_control_with_neurons"
grouping <- "celltypeloc"

plot_prefix <- paste0(dataset, "_", grouping, "_")

ppi_project <- FALSE

#######################################################################################################################################

## 1) Load the full annotated, normalized dataset  [ONLY NEED TO RUN THIS ONCE TO PREP THE DATASET]

seurat_data <- SeuratDisk::LoadH5Seurat(paste0(data_dir, dataset, ".h5Seurat"))  

# Re-level celltypeloc
seurat_data$celltypeloc <- factor(seurat_data$celltypeloc, 
                                  levels = c("neuron", "Mye", "NK", "B", "T"))

# CellChat likes to have this metadata
seurat_data$samples <- factor(paste0(seurat_data$location, "_", seurat_data$date, "_", seurat_data$condition))

cellchat <- createCellChat(object = seurat_data, group.by = "celltypeloc", assay = "SCT")

cellchat@DB <- CellChatDB.mouse # use CellChatDB.human if running on human data

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# if (ppi_project) { 
#   # project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
#   cellchat <- projectData(cellchat, PPI.mouse)
#   
#   cellchat <- computeCommunProb(cellchat, type = "triMean", raw.use = FALSE)
#   
# } else {
#   cellchat <- computeCommunProb(cellchat, type = "triMean") 
# }

#cellchat <- computeCommunProb(cellchat, type = "triMean") 
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1) 

# Users can filter out the cell-cell communication if there are only few cells in certain cell groups. 
# By default, the minimum number of cells required in each cell group for cell-cell communication is 10.
cellchat <- filterCommunication(cellchat, min.cells = 10)

# Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways

cellchat %>% readr::write_rds(paste0(data_dir, "cellchat/cellchat_", dataset, "_", grouping, ".rds"))

#######################################################################################################################################
#cellchat <- readr::read_rds(paste0(data_dir, "cellchat/cellchat_", dataset, "_", grouping, ".rds"))
cellchat <- readr::read_rds(paste0(data_dir, "cellchat/cellchat_trim10_", dataset, "_", grouping, ".rds"))

neuron_as_sender <- FALSE

## Set Source and Target Groups of Interest
if (TRUE == neuron_as_sender) {
  
  # Interested in Neuron as source
  source_groups <- match(c("neuron"), levels(cellchat@idents))
  target_groups <- match(setdiff(levels(cellchat@idents), c("neuron")), 
                         levels(cellchat@idents))
  
} else {
  
  # Interested in Neuron as receiver
  source_groups <- match(setdiff(levels(cellchat@idents), c("neuron")), 
                         levels(cellchat@idents))
  target_groups <- match(c("neuron"), levels(cellchat@idents))
}

# neurons <- match(c("neuron"), levels(cellchat@idents))
# all_other_cells <- match(c("cornea_B", "cornea_Epi", "cornea_Mye", "cornea_NK", "cornea_other", "cornea_other_imm", "cornea_T", 
#                            "tg_B", "tg_Mye", "tg_NK", "tg_other_imm", "tg_T"), 
#                          levels(cellchat@idents))

circle_plot <- 
  netVisual_circle(cellchat@net$weight, 
                   idents.use = levels(cellchat@idents),
                   targets.use = target_groups, 
                   sources.use = source_groups,
                   #vertex.weight = 2.0, 
                   vertex.label.cex = 2, 
                   #vertex.size.max = 0.1, 
                   #color.use = "darkgrey", 
                   layout = in_circle(),
                   alpha.edge = 1,
                   weight.scale = T, 
                   label.edge= T, 
                   edge.label.cex = 1.7, 
                   shape = "none",
                   edge.curved = 0.1,
                   arrow.width = 10, 
                   arrow.size = 0.2, 
                   margin = 0,
                   title.name = "Interaction weights/strength")
  
svg(paste0(figures_dir, plot_prefix, "circle_plot.svg"), width = 5, height = 5)
circle_plot
dev.off()

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
#netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, 
                 sources.use = neurons_and_tg_immune, 
                 vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

netVisual_circle(cellchat@net$weight, 
                 targets.use = neurons_and_tg_immune,
                 vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

mat <- cellchat@net$weight["neuron",]  # Just Neuron
mat <- cellchat@net$weight  # All of them
par(mfcol = c(2,3), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, 
                   vertex.weight = groupSize, 
                   weight.scale = T, 
                   edge.weight.max = max(mat), 
                   title.name = rownames(mat)[i])
}

# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")

# Heatmap
#par(mfrow=c(1,2))
cc1 <- netVisual_heatmap(cellchat, color.heatmap = "Reds", measure = "count", font.size = 18) 
cc2 <- netVisual_heatmap(cellchat, color.heatmap = "Reds", measure = "weight", font.size = 18)
cc1 + cc2


## PATHWAYS HEATMAP
source <- levels(cellchat@idents)[source_groups]
target <- levels(cellchat@idents)[target_groups]

# Show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
# signaling_poi <- cellchat@netP$pathways ## All pathways
if (TRUE == neuron_as_sender) {
  signaling_poi <- c("GALECTIN", "MIF", "CCL", "APP", "ADGRE", "PTN")
} else {
  signaling_poi <- c("GALECTIN", "APP", "CypA", "CCL", "LAMININ", "PARs", "CD39")
}

df <- data.frame()

for (pathway in cellchat@netP$pathways) {
  row <- cellchat@netP$prob[, , pathway][source, target] %>% 
    t() %>% 
    as.data.frame() %>% 
    mutate(pathway = pathway)
  
  df <- bind_rows(df, row)
}

pathway_relative_strength <- df %>% 
  column_to_rownames("pathway")

colnames(pathway_relative_strength) <- paste0(source, " -> ", target)

# Remove rows for pathways with all zeroes
pathway_relative_strength <- pathway_relative_strength %>% filter_all(any_vars(. != 0))

pathway_relative_strength <- pathway_relative_strength / max(pathway_relative_strength)

pathway_relative_strength.scaled <- scale(pathway_relative_strength)

pathway_relative_strength.log <- log1p(pathway_relative_strength)


heatmap_mtx <- pathway_relative_strength.scaled[intersect(signaling_poi, rownames(pathway_relative_strength.scaled)), ] %>% as.matrix()
pathway_heatmap.scaled <- ComplexHeatmap::Heatmap(heatmap_mtx, 
                                                  col = rev(RColorBrewer::brewer.pal(11,"Spectral")), 
                                                  cluster_rows = FALSE, 
                                                  cluster_columns = FALSE,
                                                  row_names_side = "left", 
                                                  row_title = "Pathways", row_title_gp = gpar(fontsize = 18, fontface = "bold"),
                                                  row_names_gp = gpar(fontsize = 14, fontface = "bold"),
                                                  column_names_gp = gpar(fontsize = 14, fontface = "bold"), 
                                                  column_names_rot = 45, 
                                                  heatmap_legend_param = list(title = "Commun\nProb",
                                                                              at = c(min(heatmap_mtx), max(heatmap_mtx)),
                                                                              border = "black", just = "middle",
                                                                              labels = c("min", "max")))

# pathway_heatmap.log <- 
#   ComplexHeatmap::Heatmap(pathway_relative_strength.log[signaling_poi, ] %>% as.matrix(), 
#                           col = rev(RColorBrewer::brewer.pal(11,"Spectral")), 
#                           cluster_rows = FALSE, 
#                           cluster_columns = FALSE,
#                           row_names_side = "left", 
#                           row_title = "Pathways", row_title_gp = gpar(fontsize = 18, fontface = "bold"),
#                           row_names_gp = gpar(fontsize = 14, fontface = "bold"),
#                           column_names_gp = gpar(fontsize = 14, fontface = "bold"), 
#                           column_names_rot = 45, 
#                           # heatmap_legend_param = list(title = "Commun.\nProb", 
#                           #                             #at = c(0, 1), 
#                           #                             border = "black", just = "middle",
#                           #                             labels = c("min", "max"))
#   )

svg(paste0(figures_dir, plot_prefix, "pathway_heatmap.svg"), width = 5, height = 5)
pathway_heatmap.scaled
dev.off()


## LIGAND-RECEPTOR BUBBLE PLOTS
lr_pairs <- 
  netVisual_bubble(cellchat, 
                 sources.use = source_groups, 
                 targets.use = target_groups, 
                 color.heatmap = c("Spectral"), 
                 sort.by.source.priority = T,
                 thresh = 0.05, 
                 font.size = 14, 
                 dot.size.min = 7, 
                 return.data = TRUE)$communication

sigLRs <- cellchat@LR$LRsig %>% filter(pathway_name %in% signaling_poi) %>% pull(interaction_name_2)

lr_pairs.scaled <- lr_pairs %>% 
  mutate(label = paste0(source, " -> ", target), 
         lr_pair = interaction_name_2) %>% 
  select(label, prob, lr_pair) %>% 
  pivot_wider(names_from = "label", values_from = "prob") %>% 
  mutate(across(!starts_with("lr_pair"), scale)) %>% 
  filter(lr_pair %in% sigLRs) %>%   # NOTE: Filter AFTER scaling
  pivot_longer(!lr_pair, names_to = "cc", values_to = "comm_prob") %>% 
  drop_na(comm_prob)

lr_pairs.scaled %>% 
  ggplot(., aes(x = lr_pair, y = cc, color = comm_prob)) +
  geom_point(pch = 16, size = 7) +
  theme_linedraw() + 
  guides(color = guide_colourbar(title = "Commun\nProb", ticks.colour = NA, frame.colour = "black")) + 
  scale_colour_gradientn(colors = rev(colorRampPalette(RColorBrewer::brewer.pal(n = 10, name = "Spectral"))(99)), 
                         na.value = "white", 
                         #limits=c(quantile(df$prob, 0,na.rm= T), quantile(df$prob, 1,na.rm= T)),
                         breaks = c(quantile(lr_pairs.scaled$comm_prob, 0,na.rm= T), quantile(lr_pairs.scaled$comm_prob, 1,na.rm= T)), 
                         labels = c("min","max")) +
  #guides(color = guide_colourbar(barwidth = 0.5, title = "Commun. Prob."))
  theme(legend.title = element_text(face = "bold", size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.text = element_text(face = "bold", size = 14),
        axis.title = element_blank(), 
        #aspect.ratio = 0.4, 
        panel.grid.major = element_blank()) + 
  geom_vline(xintercept=seq(1.5, length(unique(lr_pairs.scaled$lr_pair)) -0.5, 1),lwd=0.1,colour="grey90") + 
  geom_hline(yintercept=seq(1.5, length(unique(lr_pairs.scaled$cc)) -0.5, 1),lwd=0.1,colour="grey90")

ggsave(paste0(figures_dir, plot_prefix, "bubbles_select_pathways.svg"), plot = last_plot(), width = 14.5, height = 3.4)
#ggsave(paste0(figures_dir, plot_prefix, "bubbles_select_pathways_neuron_target.png"), plot = last_plot(), width = 25, height = 3.4)


top_interactions <- 
  netVisual_bubble(cellchat, 
                   sources.use = neurons, 
                   targets.use = immune, 
                   signaling = signaling_poi,
                   color.heatmap = c("Spectral"), 
                   sort.by.source.priority = T,
                   thresh = 0.05, 
                   font.size = 14, 
                   dot.size.min = 7, 
                   return.data = TRUE)
  

## Look at Ccl2 - Ccr2
#top_interactions$communication %>% filter("Ccl2" == ligand)
#cellchat.trim10@net$pval[, , "CCL2_CCR2"]

top_interactions <- 
  netVisual_bubble(cellchat, 
                   sources.use = neurons, 
                   targets.use = immune, 
                   thresh = 0.05, 
                   angle.x = 45, 
                   font.size = 14, 
                   dot.size.min = 7,
                   return.data = TRUE)

top_interactions$communication %>% 
  filter("neuron" == source) %>% 
  arrange(desc(prob)) %>% 
  slice_head(n = 25)

top_interactions$communication %>% 
  #filter("neuron" == source | "neuron" == target) %>% 
  filter("Ccl2" == ligand | "Ccl2" == receptor)


# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
png(paste0(figures_dir, plot_prefix, "_kos_re_compare_heatmap.png"), width = 20, height = 25, units = "in", res = 150)
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", width = 10, height = 25)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", width = 10, height = 25)
ht1 + ht2
dev.off()



#######################################################################################################################################

interactions <- subsetCommunication(cellchat)

# Find the Top N receptors on the neurons (allow neurons to also be senders)
interactions %>% 
  filter("neuron" == target) %>% 
  group_by(receptor) %>% 
  arrange(pval, desc(prob)) %>% 
  slice_head(n = 1) %>% 
  ungroup() %>% 
  arrange(pval, desc(prob))

# Find the Top N receptors on the neurons (do not allow neurons to also be senders)
topN_receptor_excl_neuron_send <- 
  interactions %>% 
  filter("neuron" == target, "neuron" != source) %>% 
  group_by(receptor) %>% 
  arrange(pval, desc(prob)) %>% 
  slice_head(n = 1) %>% 
  ungroup() %>% 
  arrange(pval, desc(prob))

VlnPlot(seurat_data, assay = "SCT", layer = "data", stack = TRUE, flip = TRUE,
        group.by = "celltypeloc",
        features = sub("ITGA3_\\w*", "Itga3", topN_receptor_excl_neuron_send$receptor[1:10]),
        ) + NoLegend()

# Find the Top N ligands on the neurons (allow neurons to also be receivers)
interactions %>% 
  filter("neuron" == source) %>% 
  group_by(ligand) %>% 
  arrange(pval, desc(prob)) %>% 
  slice_head(n = 1) %>% 
  ungroup() %>% 
  arrange(pval, desc(prob))

# Find the Top N ligands on the neurons (do not allow neurons to also be receivers)
topN_ligand_excl_neuron_rec <- 
  interactions %>% 
  filter("neuron" == source, "neuron" != target) %>% 
  group_by(ligand) %>% 
  arrange(pval, desc(prob)) %>% 
  slice_head(n = 1) %>% 
  ungroup() %>% 
  arrange(pval, desc(prob))

VlnPlot(seurat_data, assay = "SCT", layer = "data", stack = TRUE, flip = TRUE,
        group.by = "celltypeloc",
        features = sub("ITGA3_\\w*", "Itga3", topN_ligand_excl_neuron_rec$ligand[1:10]),
) + NoLegend()

brian_genes <- c("Ccr2", 
                 "Ifngr1", "Ifngr2", "Ifnar1", "Ifnar2", 
                 "Tnfrsf21", "Tnfrsf11a", "Tnfrsf1a", 
                 "Il4ra", "Il6", "Il6st", "Il10rb", "Il13ra1")

VlnPlot(seurat_data, assay = "SCT", layer = "data", stack = TRUE, flip = TRUE,
        group.by = "celltypeloc",
        features = brian_genes,
) + NoLegend()

#######################################################################################################################################

library(leiden)

## Neuron Markers and Co-expression with ligands from above
seurat_data.neurons <- subset(seurat_data, subset = cell_L2 == "neuron")

seurat_data.neurons <- seurat_data.neurons %>% 
  SCTransform(vst.flavor = "v2", return.only.var.genes = FALSE) %>% 
  RunPCA(assay = "SCT", npcs = 50, verbose = TRUE) %>% 
  RunUMAP(assay = "SCT", dims = 1:30) %>% 
  FindNeighbors(assay = "SCT", dims = 1:30, verbose = TRUE) %>%
  FindClusters(resolution = 1.0, method = "igraph", algorithm = "Leiden") %>% 
  PrepSCTFindMarkers(verbose = TRUE)

VlnPlot(seurat_data.neurons, 
        assay = "SCT", 
        layer = "data", 
        stack = TRUE, flip = TRUE,
        group.by = "seurat_clusters",
        features = c("Trpv1", "Calca", "Tac1", "GFRa2", "Piezo2", "Trpm8", 
                     "App", "Lgals9", "Mif", "Ccl2")) + 
  NoLegend()

FeaturePlot(seurat_data.neurons, features = c("Trpv1", "Calca", "Tac1", "GFRa2", "Piezo2", "Trpm8", 
                                              "App", "Lgals9", "Mif"))

DotPlot(seurat_data, 
        assay = "SCT", 
        group.by = "location",
        features = c("Ccl2", "Ccr2", "App"))

AverageExpression(seurat_data, 
                  assays = c("SCT"),
                  features = c("Ccl2", "Ccr2"), 
                  group.by = "location")

hist(seurat_data.neurons@assays$RNA$counts["Ccl2", ])
hist(seurat_data.neurons@assays$SCT$data["Ccl2", ])

#######################################################################################################################################
## Generate TG Immune Umap with Markers
combined_tg <- seurat_data <- SeuratDisk::LoadH5Seurat(paste0(data_dir, "combined_tg_annot.h5Seurat"))  

# Immune Cell Proportions
combined_tg_props <- combined_tg@meta.data %>% 
  filter("control" == condition) %>% 
  group_by(condition, cell_L2) %>% 
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(prop = round(n / sum(n) * 100, digits = 1))

#combined_tg$cell_L2 <- 
combined_tg$umap_label <- factor(paste0(combined_tg$cell_L2, "\n",
                                        combined_tg_props$prop[as.numeric(combined_tg$cell_L2)], 
                                        "%"))

tg_umap <- DimPlot(combined_tg, 
                   group.by = "umap_label", 
                   reduction = "umap", 
                   label = TRUE, label.box = TRUE, label.size = 5, repel = TRUE) + 
  NoLegend() + ggtitle("TG immune cell proportions")

ggsave(paste0(figures_dir, "umap.svg"), plot = tg_umap, width = 5, height = 5)

# combined_tg_props %>% 
#   ggplot(aes(cell_L2, prop)) + 
#   scale_fill_viridis_d() + 
#   geom_col(position=position_dodge()) + 
#   ggprism::theme_prism() + 
#   labs(x = "", y = "Proportion of TG immune cells") + 
#   theme(aspect.ratio = 1.7)

tg_immune_markers <- c("Ptprc", "Pax5", "Cd19", "Ighm", 
                       "Cd3d", 
                       "Klrb1c", 
                       "Adgre1", "Cd83","Cd68")

#library(patchwork)
feature_plots <- list()
for (i in tg_immune_markers) {
  feature_plots[[i]] <- 
    FeaturePlot(combined_tg, 
                features = i,
                #alpha = 0.7,
                #cols = c("darkgreen"), 
                order = TRUE, 
                max.cutoff = 1.0, 
                pt.size = 1.4) + NoLegend() + 
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
}

combined_feature_plots <- cowplot::plot_grid(plotlist = feature_plots, ncol = 5, axis = "bltr")
ggsave(paste0(figures_dir, "features.svg"), plot = combined_feature_plots, width = 10, height = 5)

#######################################################################################################################################
## Read in Peter's UMAP embeddings and use to generate Neuron UMAPs

# seurat_data <- SeuratDisk::LoadH5Seurat(paste0(data_dir, "seurat/tg_control_with_neurons.h5Seurat"))  
# 
# neurons <- subset(seurat_data, celltypeloc == "neuron")

Convert("../data_objects/Neurons_SMARTseq_Analysis.h5ad", assay = "RNA", dest = "h5seurat", overwrite = TRUE)
neurons <- LoadH5Seurat("../data_objects/Neurons_SMARTseq_Analysis.h5seurat")

