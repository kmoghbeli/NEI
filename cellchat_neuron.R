library(tidyverse)
library(ggprism)
library(Seurat)
library(CellChat)
library(patchwork)
library(ComplexHeatmap)
options(stringsAsFactors = FALSE)
conflicted::conflict_prefer_all("dplyr", quiet = TRUE)

kaveh_colors1 <- c("#a6611a", "#dfc27d", "#f5f5f5", "#80cdc1", "#018571")
kaveh_colors2 <- c("#40004b", "#762a83", "#9970ab", "#c2a5cf", "#e7d4e8", "#f7f7f7", "#d9f0d3", "#a6dba0", "#5aae61", "#1b7837", "#00441b")
kaveh_colors3 <- c("#ff0000", "#f5f5f5","#0000ff")


config_dir <- "./config/"
data_dir <- "./data/"
model_dir <- "./models/"
results_dir <- "./results/"
figures_dir <- "./figures/"

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

# source("OmicsToolbox/marker_genes.R")
# 
# seurat.tg <- subset(seurat_data, subset = location == "tg")
# other_imm_markers <- FindMarkers(seurat.tg, assay = "SCT", group.by = "cell_L2", ident.1 = "other_imm") %>% rownames_to_column("gene")
# 
# VlnPlot(seurat.tg, 
#         assay = "SCT", layer = "data", stack = TRUE, flip = TRUE,
#         group.by = "celltypeloc",
#         features = c("Ptprc", "Itgam", "Adgre1", "Itgax", marker_genes_tg),
# ) + NoLegend()
# 
# VlnPlot(seurat.tg, 
#         assay = "SCT", layer = "data", stack = TRUE, flip = TRUE,
#         group.by = "celltypeloc",
#         features = other_imm_markers$gene[1:25],
# ) + NoLegend()
#######################################################################################################################################
# neurons <- subset(seurat_data, subset = location == "corneal_afferents")
# 
# VlnPlot(neurons, features = c("Ccl2", "Ccl21a"), assay = "SCT", layer = "data")
# VlnPlot(neurons, features = c("Ccl2", "Ccl21a"), assay = "RNA", layer = "counts")
neurons <- FindVariableFeatures(neurons, assay = "SCT")
match("Ccl2", (VariableFeatures(neurons)))
match("Ccl21a", (VariableFeatures(neurons)))
match("Ccl26", (VariableFeatures(neurons)))
#######################################################################################################################################

## 1) Load the full annotated, normalized dataset

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

if (ppi_project) { 
  # project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
  cellchat <- projectData(cellchat, PPI.mouse)
  
  cellchat <- computeCommunProb(cellchat, type = "triMean", raw.use = FALSE)
  
} else {
  cellchat <- computeCommunProb(cellchat, type = "triMean") 
}

# Users can filter out the cell-cell communication if there are only few cells in certain cell groups. 
# By default, the minimum number of cells required in each cell group for cell-cell communication is 10.
cellchat <- filterCommunication(cellchat, min.cells = 10)

# Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways

cellchat %>% readr::write_rds(paste0(data_dir, "cellchat/cellchat_", dataset, "_", grouping, ".rds"))

#######################################################################################################################################

## Chat with Shams
control_tg <- subset(seurat_data, subset = location == "tg")

tg_immune_markers <- c("Ptprc", "Pax5", "Cd19", "Ighm", 
                       "Cd3d", 
                       "Klrb1c", 
                       "Adgre1", "Cd83","Cd68")

VlnPlot(control_tg, assay = "SCT", layer = "data", 
        stack = TRUE, flip = TRUE,
        group.by = "cell_L4",
        features = c(tg_immune_markers, "Ms4a1", "Cd74", "Cd79a", "Cd79b")) + NoLegend()

VlnPlot(combined_tg, assay = "SCT", layer = "data", 
        stack = TRUE, flip = TRUE,
        group.by = "condition",
        features = c("Ccr2", "Cx3cr1", "Cd69", "Itgae")) + NoLegend()

combined_tg.mye = subset(combined_tg, subset = cell_L2 == "Mye")

VlnPlot(combined_tg.mye, assay = "RNA", layer = "data", 
        stack = TRUE, flip = TRUE,
        group.by = "condition",
        features = c("Ccr2", "Cx3cr1")) + NoLegend()

combined_tg.T = subset(combined_tg, subset = cell_L2 == "T")

VlnPlot(combined_tg.T, assay = "SCT", layer = "data", 
        stack = TRUE, flip = TRUE,
        group.by = "condition",
        features = c("Cd69", "Itgae")) + NoLegend()

DimPlot(seurat_data)
          

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



#######################################################################################################################################
cellchat <- readr::read_rds(paste0(data_dir, "cellchat/cellchat_", dataset, "_", grouping, ".rds"))

# neurons <- match(c("neuron"), levels(cellchat@idents))
# all_other_cells <- match(c("cornea_B", "cornea_Epi", "cornea_Mye", "cornea_NK", "cornea_other", "cornea_other_imm", "cornea_T", 
#                            "tg_B", "tg_Mye", "tg_NK", "tg_other_imm", "tg_T"), 
#                          levels(cellchat@idents))

neurons <- match(c("neuron"), levels(cellchat@idents))
immune <- match(c("B", "Mye", "NK", "T"), levels(cellchat@idents))
neurons_and_immune <- match(c("neuron", "B", "Mye", "NK", "T"), levels(cellchat@idents))
neurons_and_nk_mye <- match(c("neuron", "Mye", "NK"), levels(cellchat@idents))

#cellchat@net$weight[c("neuron"), c("B", "Mye", "NK", "T")] %>% 
  
neuron_connect_plot <- 
  netVisual_circle(cellchat@net$weight, 
                 idents.use = levels(cellchat@idents),
                 targets.use = immune, 
                 sources.use = neurons,
                 #vertex.weight = 2.0, 
                 vertex.label.cex = 2, 
                 #vertex.size.max = 0.1, 
                 color.use = "darkgrey", 
                 layout = in_circle(),
                 alpha.edge = 1,
                 weight.scale = T, 
                 label.edge= T, edge.label.cex = 1.7, 
                 shape = "none",
                 edge.curved = 0.1,
                 arrow.width = 10, arrow.size = 0.2, margin = 0,
                 title.name = "Interaction weights/strength")

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

# Show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_bubble(cellchat, sources.use = neurons, targets.use = all_other_cells)
netVisual_bubble(cellchat, sources.use = all_other_cells, targets.use = neurons)

neurons_and_mye <- match(c("neuron", "Mye"), 
                         levels(cellchat@idents))

mye_and_nk <- match(c("Mye", "NK"), 
                    levels(cellchat@idents))

bubble_plot <- netVisual_bubble(cellchat, 
                 sources.use = neurons, 
                 targets.use = immune, 
                 color.heatmap = c("Spectral"), 
                 sort.by.source.priority = T,
                 #thresh = 0.05, 
                 angle.x = 45, 
                 font.size = 14, 
                 dot.size.min = 7)

top_interactions <- 
  netVisual_bubble(cellchat, 
                   sources.use = neurons_and_mye, 
                   targets.use = mye_and_nk, 
                   thresh = 0.05, 
                   angle.x = 45, 
                   font.size = 14, 
                   dot.size.min = 7,
                   return.data = TRUE
  )

top_interactions$communication %>% 
  filter("neuron" == source | "neuron" == target) %>% 
  arrange(desc(prob)) %>% 
  slice_head(n = 25)

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
png(paste0(figures_dir, plot_prefix, "_kos_re_compare_heatmap.png"), width = 20, height = 25, units = "in", res = 150)
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", width = 10, height = 25)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", width = 10, height = 25)
ht1 + ht2
dev.off()

df <- data.frame()

for (pathway in cellchat@netP$pathways) {
  row <- cellchat@netP$prob[, , pathway]["neuron", immune] %>% 
    t() %>% 
    as.data.frame() %>% 
    mutate(pathway = pathway)

  df <- bind_rows(df, row)
}

pathway_relative_strength <- df %>% 
  column_to_rownames("pathway")

pathway_relative_strength <- pathway_relative_strength / max(pathway_relative_strength)
colnames(pathway_relative_strength) <- paste0("Neuron -> ", colnames(pathway_relative_strength))

# Remove rows for pathways with all zeroes
pathway_relative_strength <- pathway_relative_strength %>% filter_all(any_vars(. != 0))

pathway_heatmap <- 
  ComplexHeatmap::Heatmap(pathway_relative_strength, 
                          col = rev(RColorBrewer::brewer.pal(11,"Spectral")), 
                          cluster_rows = FALSE, 
                          cluster_columns = FALSE,
                          row_names_side = "left", 
                          row_title = "Pathways", row_title_gp = gpar(fontsize = 18, fontface = "bold"),
                          row_names_gp = gpar(fontsize = 14, fontface = "bold"),
                          column_names_gp = gpar(fontsize = 14, fontface = "bold"), 
                          column_names_rot = 45, 
                          heatmap_legend_param = list(title = "Commun.\nProb", at = c(0, 1), 
                                                      border = "black", just = "middle",
                                                      labels = c("min", "max"))
  )

## Save Plots
ggsave(paste0(figures_dir, "umap.svg"), plot = tg_umap, width = 5, height = 5)
ggsave(paste0(figures_dir, "features.svg"), plot = combined_feature_plots, width = 10, height = 5)
svg(paste0(figures_dir, "neuron_circle.svg"), width = 5, height = 5)
neuron_connect_plot
dev.off()
svg(paste0(figures_dir, "pathway_heatmap.svg"), width = 5, height = 5)
pathway_heatmap
dev.off()
ggsave(paste0(figures_dir, "bubbles.svg"), plot = bubble_plot, width = 5, height = 7)


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
                     "App", "Lgals9", "Mif")
) + NoLegend()

FeaturePlot(seurat_data.neurons, features = c("Trpv1", "Calca", "Tac1", "GFRa2", "Piezo2", "Trpm8", 
                                              "App", "Lgals9", "Mif"))

#######################################################################################################################################

netVisual_circleKM <-function(net, color.use = NULL,title.name = NULL, sources.use = NULL, targets.use = NULL, idents.use = NULL, remove.isolate = FALSE, top = 1,
                            weight.scale = FALSE, vertex.weight = 20, vertex.weight.max = NULL, vertex.size.max = NULL, vertex.label.cex=1,vertex.label.color= "black",
                            edge.weight.max = NULL, edge.width.max=8, alpha.edge = 0.6, label.edge = FALSE,edge.label.color='black',edge.label.cex=0.8,
                            edge.curved=0.2,shape='circle',layout=in_circle(), margin=0.2, vertex.size = NULL,
                            arrow.width=1,arrow.size = 0.2,
                            text.x = 0, text.y = 1.5){
  if (!is.null(vertex.size)) {
    warning("'vertex.size' is deprecated. Use `vertex.weight`")
  }
  if (is.null(vertex.size.max)) {
    if (length(unique(vertex.weight)) == 1) {
      vertex.size.max <- 5
    } else {
      vertex.size.max <- 15
    }
  }
  options(warn = -1)
  thresh <- stats::quantile(net, probs = 1-top)
  net[net < thresh] <- 0
  
  if ((!is.null(sources.use)) | (!is.null(targets.use)) | (!is.null(idents.use)) ) {
    if (is.null(rownames(net))) {
      stop("The input weighted matrix should have rownames!")
    }
    cells.level <- rownames(net)
    df.net <- reshape2::melt(net, value.name = "value")
    colnames(df.net)[1:2] <- c("source","target")
    # keep the interactions associated with sources and targets of interest
    if (!is.null(sources.use)){
      if (is.numeric(sources.use)) {
        sources.use <- cells.level[sources.use]
      }
      df.net <- subset(df.net, source %in% sources.use)
    }
    if (!is.null(targets.use)){
      if (is.numeric(targets.use)) {
        targets.use <- cells.level[targets.use]
      }
      df.net <- subset(df.net, target %in% targets.use)
    }
    if (!is.null(idents.use)) {
      if (is.numeric(idents.use)) {
        idents.use <- cells.level[idents.use]
      }
      df.net <- filter(df.net, (source %in% idents.use) | (target %in% idents.use))
    }
    df.net$source <- factor(df.net$source, levels = cells.level)
    df.net$target <- factor(df.net$target, levels = cells.level)
    df.net$value[is.na(df.net$value)] <- 0
    net <- tapply(df.net[["value"]], list(df.net[["source"]], df.net[["target"]]), sum)
  }
  net[is.na(net)] <- 0
  
  if (is.null(color.use)) {
    color.use = scPalette(nrow(net))
    names(color.use) <- rownames(net)
  } else {
    if (is.null(names(color.use))) {
      stop("The input `color.use` should be a named vector! \n")
    }
  }
  if (remove.isolate) {
    idx1 <- which(Matrix::rowSums(net) == 0)
    idx2 <- which(Matrix::colSums(net) == 0)
    idx.isolate <- intersect(idx1, idx2)
    if (length(idx.isolate) > 0) {
      net <- net[-idx.isolate, ]
      net <- net[, -idx.isolate]
      color.use = color.use[-idx.isolate]
      if (length(unique(vertex.weight)) > 1) {
        vertex.weight <- vertex.weight[-idx.isolate]
      }
    }
  }
  
  g <- graph_from_adjacency_matrix(net, mode = "directed", weighted = T)
  edge.start <- igraph::ends(g, es=igraph::E(g), names=FALSE)
  coords<-layout_(g,layout)
  if(nrow(coords)!=1){
    coords_scale=scale(coords)
  }else{
    coords_scale<-coords
  }
  
  if (is.null(vertex.weight.max)) {
    vertex.weight.max <- max(vertex.weight)
  }
  #vertex.weight <- vertex.weight/vertex.weight.max*vertex.size.max+5
  vertex.weight <- 10
  
  loop.angle<-ifelse(coords_scale[igraph::V(g),1]>0,-atan(coords_scale[igraph::V(g),2]/coords_scale[igraph::V(g),1]),pi-atan(coords_scale[igraph::V(g),2]/coords_scale[igraph::V(g),1]))
  igraph::V(g)$size<-vertex.weight
  igraph::V(g)$color<-color.use[igraph::V(g)]
  igraph::V(g)$frame.color <- color.use[igraph::V(g)]
  igraph::V(g)$label.color <- vertex.label.color
  igraph::V(g)$label.cex<-vertex.label.cex
  if(label.edge){
    igraph::E(g)$label<-igraph::E(g)$weight
    igraph::E(g)$label <- round(igraph::E(g)$label, digits = 1)
  }
  if (is.null(edge.weight.max)) {
    edge.weight.max <- max(igraph::E(g)$weight)
  }
  
  # if (weight.scale == TRUE) {
  #   #E(g)$width<-0.3+edge.width.max/(max(E(g)$weight)-min(E(g)$weight))*(E(g)$weight-min(E(g)$weight))
  #   igraph::E(g)$width<- 0.3+igraph::E(g)$weight/edge.weight.max*edge.width.max
  # }else{
  #   igraph::E(g)$width<-0.3+edge.width.max*igraph::E(g)$weight
  # }
  
  igraph::E(g)$width <- 1
  
  igraph::E(g)$arrow.width<-arrow.width
  igraph::E(g)$arrow.size<-arrow.size
  igraph::E(g)$label.color<-edge.label.color
  igraph::E(g)$label.cex<-edge.label.cex
  #igraph::E(g)$color<- grDevices::adjustcolor(igraph::V(g)$color[edge.start[,1]],alpha.edge)
  igraph::E(g)$color<- "blue"
  igraph::E(g)$loop.angle <- rep(0, length(igraph::E(g)))
  
  if(sum(edge.start[,2]==edge.start[,1])!=0){
    igraph::E(g)$loop.angle[which(edge.start[,2]==edge.start[,1])]<-loop.angle[edge.start[which(edge.start[,2]==edge.start[,1]),1]]
  }
  radian.rescale <- function(x, start=0, direction=1) {
    c.rotate <- function(x) (x + start) %% (2 * pi) * direction
    c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
  }
  label.locs <- radian.rescale(x=1:length(igraph::V(g)), direction=-1, start=0)
  label.dist <- vertex.weight/max(vertex.weight)+2
  plot(g,edge.curved=edge.curved,vertex.shape=shape,layout=coords_scale,margin=margin, vertex.label.dist=label.dist,
       vertex.label.degree=label.locs, vertex.label.family="Helvetica", edge.label.family="Helvetica") # "sans"
  if (!is.null(title.name)) {
    text(text.x,text.y,title.name, cex = 1.1)
  }
  # https://www.andrewheiss.com/blog/2016/12/08/save-base-graphics-as-pseudo-objects-in-r/
  # grid.echo()
  # gg <-  grid.grab()
  gg <- recordPlot()
  return(gg)
}

