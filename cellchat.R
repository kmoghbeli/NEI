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

#dataset <- "sept_cornea"
dataset <- "combined_tg"
grouping <- "cell_L2"
#grouping <- "cell_L3"
#grouping <- "cell_L4"

plot_prefix <- paste0(dataset, "_", grouping, "_")

## 1) Load the full annotated, normalized dataset

seurat_data <- SeuratDisk::LoadH5Seurat(paste0(data_dir, dataset, "_annot.h5Seurat"))  

# CellChat likes to have this metadata
seurat_data$samples <- factor(paste0(seurat_data$location, "_", seurat_data$date, "_", seurat_data$condition))

# Create the subsets for each condition
seurat_data.control <- subset(seurat_data, subset = "control" == condition)
seurat_data.scratch <- subset(seurat_data, subset = "scratch" == condition)
seurat_data.kos <- subset(seurat_data, subset = "kos" == condition)
seurat_data.re <- subset(seurat_data, subset = "re" == condition)

# Re-level the groupings (sometimes subsets don't contain all labels/values and then CellChat complains)
seurat_data.control <- AddMetaData(seurat_data.control, 
                                   seurat_data.control[[grouping]] %>% 
                                     rownames_to_column() %>% deframe() %>% 
                                     factor() %>% 
                                     enframe("barcode", grouping) %>% 
                                     column_to_rownames("barcode"))

seurat_data.scratch <- AddMetaData(seurat_data.scratch, 
                                   seurat_data.scratch[[grouping]] %>% 
                                     rownames_to_column() %>% deframe() %>% 
                                     factor() %>% 
                                     enframe("barcode", grouping) %>% 
                                     column_to_rownames("barcode"))

seurat_data.kos <- AddMetaData(seurat_data.kos, 
                               seurat_data.kos[[grouping]] %>% 
                                 rownames_to_column() %>% deframe() %>% 
                                 factor() %>% 
                                 enframe("barcode", grouping) %>% 
                                 column_to_rownames("barcode"))

seurat_data.re <- AddMetaData(seurat_data.re, 
                              seurat_data.re[[grouping]] %>% 
                                rownames_to_column() %>% deframe() %>% 
                                factor() %>% 
                                enframe("barcode", grouping) %>% 
                                column_to_rownames("barcode"))

## From the CellChat tutorial: 
## https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/CellChat-vignette.html

cellchat.control <- createCellChat(object = seurat_data.control, group.by = grouping, assay = "SCT")
cellchat.scratch <- createCellChat(object = seurat_data.scratch, group.by = grouping, assay = "SCT")
cellchat.kos <- createCellChat(object = seurat_data.kos, group.by = grouping, assay = "SCT")
cellchat.re <- createCellChat(object = seurat_data.re, group.by = grouping, assay = "SCT")

CellChatDB <- CellChatDB.mouse # use CellChatDB.human if running on human data
showDatabaseCategory(CellChatDB)

dplyr::glimpse(CellChatDB$interaction)

# use a subset of CellChatDB for cell-cell communication analysis
# CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling
# Only uses the Secreted Signaling from CellChatDB v1
#  CellChatDB.use <- subsetDB(CellChatDB, search = list(c("Secreted Signaling"), c("CellChatDB v1")), key = c("annotation", "version"))
# use all CellChatDB except for "Non-protein Signaling" for cell-cell communication analysis
# CellChatDB.use <- subsetDB(CellChatDB)

# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB. We do not suggest to use it in this way because CellChatDB v2 includes "Non-protein Signaling" (i.e., metabolic and synaptic signaling). 

# set the used database in the object
cellchat.control@DB <- CellChatDB.use
cellchat.scratch@DB <- CellChatDB.use
cellchat.kos@DB <- CellChatDB.use
cellchat.re@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cellchat.control <- subsetData(cellchat.control) # This step is necessary even if using the whole database
cellchat.scratch <- subsetData(cellchat.scratch)
cellchat.kos <- subsetData(cellchat.kos)
cellchat.re <- subsetData(cellchat.re)

future::plan("multisession", workers = 4) # do parallel
options(future.globals.maxSize = 1000 * 1024^2)   ## 1GB
options(future.rng.onMisuse = "ignore")

cellchat.control <- identifyOverExpressedGenes(cellchat.control)
cellchat.scratch <- identifyOverExpressedGenes(cellchat.scratch)
cellchat.kos <- identifyOverExpressedGenes(cellchat.kos)
cellchat.re <- identifyOverExpressedGenes(cellchat.re)

cellchat.control <- identifyOverExpressedInteractions(cellchat.control)
cellchat.scratch <- identifyOverExpressedInteractions(cellchat.scratch)
cellchat.kos <- identifyOverExpressedInteractions(cellchat.kos)
cellchat.re <- identifyOverExpressedInteractions(cellchat.re)

# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
# cellchat <- projectData(cellchat, PPI.human)

# Compute the communication probability and infer cellular communication network
# (See: Part II: Inference of cell-cell communication network in the Tutorial for more info on tuning here)
# By default type = "triMean", producing fewer but stronger interactions. 
# When setting type = "truncatedMean", a value should be assigned to trim, producing more interactions.
# To determine a proper value of trim, CellChat provides a function computeAveExpr, which can help to check 
# the average expression of signaling genes of interest, e.g,: 
# computeAveExpr(cellchat, features = c("CXCL12","CXCR4"), type =  "truncatedMean", trim = 0.1) # (Note gene names should be appropriate for organism)

#cellchat <- computeCommunProb(cellchat, type = "triMean")  # This is equivalent to using truncatedMean with trim of 0.25 (i.e., MORE stringent)
#cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1)  

cellchat.control <- computeCommunProb(cellchat.control, type = "triMean")
cellchat.scratch <- computeCommunProb(cellchat.scratch, type = "triMean")
cellchat.kos <- computeCommunProb(cellchat.kos, type = "triMean")
cellchat.re <- computeCommunProb(cellchat.re, type = "triMean")

# Users can filter out the cell-cell communication if there are only few cells in certain cell groups. 
# By default, the minimum number of cells required in each cell group for cell-cell communication is 10.
cellchat.control <- filterCommunication(cellchat.control, min.cells = 10)
cellchat.scratch <- filterCommunication(cellchat.scratch, min.cells = 10)
cellchat.kos <- filterCommunication(cellchat.kos, min.cells = 10)
cellchat.re <- filterCommunication(cellchat.re, min.cells = 10)

# Extract the inferred cellular communication network as a data frame
# We provide a function subsetCommunication to easily access the inferred cell-cell communications of interest. For example,
# df.net <- subsetCommunication(cellchat) returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. Set slot.name = "netP" to access the the inferred communications at the level of signaling pathways
# df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5)) gives the inferred cell-cell communications sending from cell groups 1 and 2 to cell groups 4 and 5.
# df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb")) gives the inferred cell-cell communications mediated by signaling WNT and TGFb.

# Infer the cell-cell communication at a signaling pathway level
cellchat.control <- computeCommunProbPathway(cellchat.control)
cellchat.scratch <- computeCommunProbPathway(cellchat.scratch)
cellchat.kos <- computeCommunProbPathway(cellchat.kos)
cellchat.re <- computeCommunProbPathway(cellchat.re)

# Calculate the aggregated cell-cell communication network
cellchat.control <- aggregateNet(cellchat.control)
cellchat.scratch <- aggregateNet(cellchat.scratch)
cellchat.kos <- aggregateNet(cellchat.kos)
cellchat.re <- aggregateNet(cellchat.re)

# Compute the network centrality scores
cellchat.control <- netAnalysis_computeCentrality(cellchat.control)
cellchat.scratch <- netAnalysis_computeCentrality(cellchat.scratch)
cellchat.kos <- netAnalysis_computeCentrality(cellchat.kos)
cellchat.re <- netAnalysis_computeCentrality(cellchat.re)

# Save
cellchat.control %>% readr::write_rds(paste0(data_dir, "cellchat/", dataset, "_cellchat_control_", grouping, ".rds"))
cellchat.scratch %>% readr::write_rds(paste0(data_dir, "cellchat/", dataset, "_cellchat_scratch_", grouping, ".rds"))
cellchat.kos %>% readr::write_rds(paste0(data_dir, "cellchat/", dataset, "_cellchat_kos_", grouping, ".rds"))
cellchat.re %>% readr::write_rds(paste0(data_dir, "cellchat/", dataset, "_cellchat_re_", grouping, ".rds"))


# # CellChat can also visualize the aggregated cell-cell communication network. 
# # For example, showing the number of interactions or the total interaction strength (weights) between any two cell groups using circle plot.
# groupSize <- as.numeric(table(cellchat@idents))
# par(mfrow = c(1,2), xpd=TRUE)
# netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
# netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
# 
# mat <- cellchat@net$weight
# par(mfrow = c(3,4), xpd=TRUE)
# for (i in 1:nrow(mat)) {
#   mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
#   mat2[i, ] <- mat[i, ]
#   netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
# }
# 
# 
# netVisual_bubble(cellchat, 
#                  sources.use = match("T", levels(cellchat@idents)), 
#                  targets.use = c(1:length(levels(cellchat@idents))), 
#                  #signaling = c("CCL","CXCL"), 
#                  remove.isolate = TRUE)
# 
# plotGeneExpression(cellchat, signaling = "CXCL", enriched.only = TRUE, type = "violin")


#######################################################################################################################################

### Visualizations ###

dataset <- "combined_tg"
#dataset <- "sept_cornea"
grouping <- "cell_L2"
#grouping <- "cell_L3"
#grouping <- "cell_L4"

plot_prefix <- paste0(dataset, "_", grouping, "_")

## Comparison between KOS and RE

cellchat.kos <- readr::read_rds(paste0(data_dir, "cellchat/", dataset, "_cellchat_kos_", grouping, ".rds"))
cellchat.re <- readr::read_rds(paste0(data_dir, "cellchat/", dataset, "_cellchat_re_", grouping, ".rds"))

## Set Source and Target Groups of Interest
if (grepl("tg", dataset, ignore.case = TRUE)) {
  
  # Tg: Mye -> Mye
  
  source_groups <- match(c("Mye"), levels(cellchat.kos@idents))
  target_groups <- match(c("Mye"), levels(cellchat.kos@idents))
  
} else if(grepl("cornea", dataset, ignore.case = TRUE)) {
  
  # Cornea: NK -> Mye, B -> Mye, B -> Epi, T->T, T->Epi
  
  source_groups <- match(c("NK", "B", "T"), levels(cellchat.kos@idents))
  target_groups <- match(c("Mye", "T", "Epi"), levels(cellchat.kos@idents))
}

object.list <- list(kos = cellchat.kos, re = cellchat.re)

cellchat <- mergeCellChat(object.list, add.names = names(object.list))

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2)) + GG_KM_THEME + theme(legend.position="none")
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight") + GG_KM_THEME + theme(legend.position="none")
ggsave(paste0(figures_dir, plot_prefix, "_kos_re_compare_int.png") , plot = gg1 + gg2, width = 5, height = 7)


# Circle plot showing differential number of interactions or interaction strength among different cell populations across two datasets
# The differential number of interactions or interaction strength in the cell-cell communication 
# network between two datasets can be visualized using circle plot, 
# where red (or blue) colored edges represent increased (or decreased) signaling in the 
# second dataset compared to the first one.
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

# Heatmap showing differential number of interactions or interaction strength among different cell populations across two datasets
# The top colored bar plot represents the sum of each column of the absolute values displayed in the 
# heatmap (incoming signaling). The right colored bar plot represents the sum of each row of the absolute values 
# (outgoing signaling). Therefore, the bar height indicates the degree of change in terms of the number of interactions
# or interaction strength between the two conditions. In the colorbar, red (or blue) represents increased (or decreased) 
# signaling in the second dataset compared to the first one.

png(paste0(figures_dir, plot_prefix, "_kos_re_compare_heatmap.png"), width = 10, height = 7, units = "in", res = 150)
gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2
dev.off()

weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}



# Identify cell populations with significant changes in sending or receiving signals
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax) + GG_KM_THEME
}
p <- cowplot::plot_grid(plotlist = gg, align = "hv", axis = "b")
ggsave(paste0(figures_dir, plot_prefix, "_kos_re_compare_cell_pops.png"), plot = p, width = 10, height = 10)


# Identify the signaling changes of specific cell populations
# gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Mac") #, signaling.exclude = "MIF")
# gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "T") #, signaling.exclude = c("MIF"))
# 
# p <- patchwork::wrap_plots(plots = list(gg1,gg2))
# ggsave(paste0(figures_dir, plot_prefix, "_kos_re_compare_Mac_T.png"), plot = p)

# Identify signaling groups based on their functional similarity
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, type = "functional")
#> Manifold learning of the signaling networks for datasets 1 2
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
#> 2D visualization of signaling networks from datasets 1 2

# Identify altered signaling with distinct interaction strength
# (A) Compare the overall information flow of each signaling pathway or ligand-receptor pair
# CellChat can identify the conserved and context-specific signaling pathways by 
# simply comparing the information flow for each signaling pathway, which is defined by 
# the sum of communication probability among all pairs of cell groups in the inferred network (i.e., the total weights in the network).
gg1 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE) + GG_KM_THEME
gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE) + GG_KM_THEME
gg1 + gg2
ggsave(paste0(figures_dir, plot_prefix, "_kos_re_inf_flow.png") , plot = gg1 + gg2, width = 14, height = 14)

# (B) Compare outgoing (or incoming) signaling patterns associated with each cell population
i = 1
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 7, height = 16)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 7, height = 16)
png(paste0(figures_dir, plot_prefix, "_kos_re_outgoing.png"), width = 10, height = 10, units = "in", res = 150)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

ht3 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 7, height = 16, color.heatmap = "GnBu")
ht4 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 7, height = 16, color.heatmap = "GnBu")
png(paste0(figures_dir, plot_prefix, "_kos_re_incoming.png"), width = 10, height = 10, units = "in", res = 300)
draw(ht3 + ht4) #, ht_gap = unit(0.5, "cm"))
dev.off()

ht5 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 7, height = 14, color.heatmap = "OrRd")
ht6 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 7, height = 14, color.heatmap = "OrRd")
png(paste0(figures_dir, plot_prefix, "_kos_re_out_and_in.png"), width = 10, height = 10, units = "in", res = 150)
draw(ht5 + ht6, ht_gap = unit(0.5, "cm"))
dev.off()

# Identify dysfunctional signaling by comparing the communication probabities
# CellChat can compare the communication probabilities mediated by L-R pairs from certain cell 
# groups to other cell groups. This can be done by setting comparison in the function netVisual_bubble.
netVisual_bubble(cellchat, sources.use = source_groups, targets.use = target_groups,  comparison = c(1, 2), angle.x = 45)
ggsave(paste0(figures_dir, plot_prefix, "_kos_re_mac_t_pathways_all.png"), width = 10, height = 7)


gg1 <- netVisual_bubble(cellchat, sources.use = source_groups, targets.use = target_groups,  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in RE", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = source_groups, targets.use = target_groups,  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in RE", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2
ggsave(paste0(figures_dir, plot_prefix, "_kos_re_mac_t_pathways_up_down.png"), width = 14, height = 7)


# Identify dysfunctional signaling by using differential expression analysis
# The above method for identifying the upgulated and down-regulated signaling is perfomed by comparing the communication probability between two datasets for each L-R pair and each pair of cell groups. Alternative, we can identify the upgulated and down-regulated signaling ligand-receptor pairs based on the differential expression analysis (DEA). Specifically, we perform differential expression analysis between two biological conditions (i.e., NL and LS) for each cell group, and then obtain the upgulated and down-regulated signaling based on the fold change of ligands in the sender cells and receptors in the receiver cells.
# Of note, users may observe the same LR pairs appearing in both the up-regulated and down-regulated results due to the fact that DEA between conditions is performed for each cell group. To perform DEA between conditions by ignoring cell group information, users can set group.DE.combined = TRUE in identifyOverExpressedGenes for CellChat v2.1.1.
# define a positive dataset, i.e., the dataset with positive fold change against the other dataset

# perform differential expression analysis 
# Of note, compared to CellChat version < v2, CellChat v2 now performs an ultra-fast Wilcoxon test using the presto package, which gives smaller values of logFC. Thus we here set a smaller value of thresh.fc compared to the original one (thresh.fc = 0.1). Users can also provide a vector and dataframe of customized DEGs by modifying the cellchat@var.features$LS.merged and cellchat@var.features$LS.merged.info. 

cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", 
                                       pos.dataset = "re",
                                       only.pos = FALSE, 
                                       thresh.pc = 0.1, 
                                       thresh.fc = 0.05, 
                                       thresh.p = 0.05, 
                                       group.DE.combined = FALSE) 

#> Use the joint cell labels from the merged CellChat object

# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = "features", variable.all = TRUE)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat, net = net, datasets = "re", ligand.logFC = 0.05, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated receptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat, net = net, datasets = "kos", ligand.logFC = -0.05, receptor.logFC = NULL)

gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

# Users can also find all the significant outgoing/incoming/both signaling according to the customized features and cell groups of interest
df <- findEnrichedSignaling(object.list[[2]], features = c("CCL19", "CXCL12"), idents = c("Inflam. FIB", "COL11A1+ FIB"), pattern ="outgoing")

# (A) Bubble plot
# We then visualize the upgulated and down-regulated signaling ligand-receptor pairs using bubble plot or chord diagram.

pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, 
                        pairLR.use = pairLR.use.up, sources.use = source_groups, targets.use = target_groups, comparison = c(1, 2),  angle.x = 90, remove.isolate = T, 
                        color.text = c("blue", "red"), 
                        title.name = paste0(str_to_title(dataset), ": Up-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = source_groups, targets.use = target_groups, comparison = c(1, 2),  angle.x = 90, remove.isolate = T, 
                        color.text = c("blue", "red"), 
                        title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
gg1 + gg2
ggsave(paste0(figures_dir, plot_prefix, "_kos_re_mac_t_pathways_up_de.png"), plot = gg1, width = 7, height = 7)
ggsave(paste0(figures_dir, plot_prefix, "_kos_re_mac_t_pathways_up_down_de.png"), width = 14, height = 7)


# pathways.show <- c("CXCL") 
# weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
# par(mfrow = c(1,2), xpd=TRUE)
# for (i in 1:length(object.list)) {
#   netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
# }
# 
# 
# cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("NL", "LS")) # set factor level
# plotGeneExpression(cellchat, signaling = "CXCL", split.by = "datasets", colors.ggplot = T, type = "violin")

#######################################################################################################################################

