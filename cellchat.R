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

# For CellChat Parallel Processing
future::plan("multisession", workers = 4) # do parallel
options(future.globals.maxSize = 1000 * 1024^2)   ## 1GB
options(future.rng.onMisuse = "ignore")

#######################################################################################################################################
#### ONLY NEED TO RUN THIS ONCE!!! ####

conditions <- c("kos", "re")
trimCutoff = 0.1
grouping <- "cell_L2"   # This is the celltype metadata variable  

# Read in Neurons
neurons <- SeuratDisk::LoadH5Seurat(paste0(data_dir, "seurat/neurons_all_conditions.h5Seurat"))

## Create TG+Neuron and Cornea+Neuron disease-specific seurat objects (to be used by CellOracle later)
for (dataset in c("combined_tg", "sept_cornea")) {
  
  message("Processing neurons + ", dataset)
  
  seurat_obj <- SeuratDisk::LoadH5Seurat(paste0(data_dir, "seurat/", dataset, "_annot.h5Seurat"))
  seurat_obj$celltypeloc <- seurat_obj$cell_L2
  
  for (cond in conditions) {
    
    message("Condition: ", cond)
    
    neurons.condition <- subset(neurons, subset = cond == condition)
    neurons.condition <- PrepSCTFindMarkers(neurons.condition)
    
    seurat_obj.condition <- subset(seurat_obj, subset = cond == condition)
    seurat_obj.condition <- PrepSCTFindMarkers(seurat_obj.condition)
    
    seurat_obj.condition <- merge(x = neurons.condition, y = c(seurat_obj.condition), merge.data = TRUE)
    
    # Re-level the groupings (sometimes subsets don't contain all labels/values and then CellChat complains)
    seurat_obj.condition <- AddMetaData(seurat_obj.condition, 
                                        seurat_obj.condition[[grouping]] %>% 
                                          rownames_to_column() %>% deframe() %>% 
                                          factor() %>% 
                                          enframe("barcode", grouping) %>% 
                                          column_to_rownames("barcode"))
    
    # CellChat likes to have this metadata
    seurat_obj.condition$samples <- factor(paste0(seurat_obj.condition$location, "_", 
                                                  seurat_obj.condition$date, "_", 
                                                  seurat_obj.condition$condition))
    
    ## https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/CellChat-vignette.html
    
    ### Creation of CellChat object
    
    cellchat.condition <- createCellChat(object = seurat_obj.condition, group.by = grouping, assay = "SCT")
    
    # showDatabaseCategory(CellChatDB)
    # 
    # dplyr::glimpse(CellChatDB$interaction)
    
    # use a subset of CellChatDB for cell-cell communication analysis
    # CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling
    # Only uses the Secreted Signaling from CellChatDB v1
    #  CellChatDB.use <- subsetDB(CellChatDB, search = list(c("Secreted Signaling"), c("CellChatDB v1")), key = c("annotation", "version"))
    # use all CellChatDB except for "Non-protein Signaling" for cell-cell communication analysis
    # CellChatDB.use <- subsetDB(CellChatDB)
    
    cellchat.condition@DB <- CellChatDB.mouse
    
    # subset the expression data of signaling genes for saving computation cost
    cellchat.condition <- subsetData(cellchat.condition) # This step is necessary even if using the whole database
    
    cellchat.condition <- identifyOverExpressedGenes(cellchat.condition)
    
    cellchat.condition <- identifyOverExpressedInteractions(cellchat.condition)
    
    # project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
    # cellchat.condition <- projectData(cellchat.condition, PPI.mouse)
    
    # Compute the communication probability and infer cellular communication network
    # (See: Part II: Inference of cell-cell communication network in the Tutorial for more info on tuning here)
    # By default type = "triMean", producing fewer but stronger interactions. 
    # When setting type = "truncatedMean", a value should be assigned to trim, producing more interactions.
    # To determine a proper value of trim, CellChat provides a function computeAveExpr, which can help to check 
    # the average expression of signaling genes of interest, e.g,: 
    # computeAveExpr(cellchat, features = c("CXCL12","CXCR4"), type =  "truncatedMean", trim = 0.1) # (Note gene names should be appropriate for organism)
    
    #cellchat.condition <- computeCommunProb(cellchat.condition, type = "triMean")  # This is equivalent to using truncatedMean with trim of 0.25 (i.e., MORE stringent)
    #cellchat.condition <- computeCommunProb(cellchat.condition, type = "truncatedMean", trim = 0.1)  
    
    cellchat.condition <- computeCommunProb(cellchat.condition, type = "truncatedMean", trim = trimCutoff)
    
    # Users can filter out the cell-cell communication if there are only few cells in certain cell groups. 
    # By default, the minimum number of cells required in each cell group for cell-cell communication is 10.
    cellchat.condition <- filterCommunication(cellchat.condition, min.cells = 10)
    
    # Infer the cell-cell communication at a signaling pathway level
    cellchat.condition <- computeCommunProbPathway(cellchat.condition)
    
    # Calculate the aggregated cell-cell communication network
    cellchat.condition <- aggregateNet(cellchat.condition) 
    
    # Compute the network centrality scores
    cellchat.condition <- netAnalysis_computeCentrality(cellchat.condition)
    
    # Save
    cellchat.condition %>% readr::write_rds(paste0(data_dir, "cellchat/cellchat_", dataset , "_neuron_", cond, "_", grouping, "_trim", trimCutoff, ".rds"))
  }
}

#######################################################################################################################################
### Visualizations ###
## Comparison between KOS and RE

#dataset <- "sept_cornea"
dataset <- "combined_tg"
grouping <- "cell_L2"
neuron_as_sender <- FALSE
#grouping <- "cell_L3"
#grouping <- "cell_L4"
trimCutoff = 0.1

cellchat.kos <- readr::read_rds(paste0(data_dir, "cellchat/cellchat_", dataset, "_neuron_kos_", grouping, "_trim", trimCutoff, ".rds"))
cellchat.re <- readr::read_rds(paste0(data_dir, "cellchat/cellchat_", dataset, "_neuron_re_", grouping, "_trim", trimCutoff, ".rds"))
object.list <- list(kos = cellchat.kos, re = cellchat.re)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

plot_prefix <- paste0("cellchat_", dataset, "_neuron_kosVre_", grouping, "_trim", trimCutoff)

## Set Source and Target Groups of Interest
if (TRUE == neuron_as_sender) {
  
  # Interested in Neuron as source
  source_groups <- match(c("neuron"), levels(cellchat.kos@idents))
  target_groups <- match(setdiff(levels(cellchat.kos@idents), c("neuron")), 
                         levels(cellchat.kos@idents))
  
} else {
  
  # Interested in Neuron as receiver
  source_groups <- match(setdiff(levels(cellchat.kos@idents), c("neuron")), 
                         levels(cellchat.kos@idents))
  target_groups <- match(c("neuron"), levels(cellchat.kos@idents))
}



# gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2)) + GG_KM_THEME + theme(legend.position="none")
# gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight") + GG_KM_THEME + theme(legend.position="none")
# gg1 + gg2
# ggsave(paste0(figures_dir, plot_prefix, "_kos_re_compare_int.png") , plot = gg1 + gg2, width = 5, height = 7)


# Circle plot showing differential number of interactions or interaction strength among different cell populations across two datasets
# The differential number of interactions or interaction strength in the cell-cell communication 
# network between two datasets can be visualized using circle plot, 
# where red (or blue) colored edges represent increased (or decreased) signaling in the 
# second dataset compared to the first one.
#png(paste0(figures_dir, plot_prefix, "_celltype_circle.png"), width = 10, height = 7, units = "in", res = 150)
par(mfrow = c(1,2), xpd = TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, sources.use = source_groups, targets.use = target_groups)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight", sources.use = source_groups, targets.use = target_groups)
#dev.off()

# Heatmap showing differential number of interactions or interaction strength among different cell populations across two datasets
# The top colored bar plot represents the sum of each column of the absolute values displayed in the 
# heatmap (incoming signaling). The right colored bar plot represents the sum of each row of the absolute values 
# (outgoing signaling). Therefore, the bar height indicates the degree of change in terms of the number of interactions
# or interaction strength between the two conditions. In the colorbar, red (or blue) represents increased (or decreased) 
# signaling in the second dataset compared to the first one.

png(paste0(figures_dir, plot_prefix, "_compare_heatmap.png"), width = 10, height = 7, units = "in", res = 150)
gg1 <- netVisual_heatmap(cellchat, comparison = c(1,2))
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2
dev.off()


# Identify cell populations with significant changes in sending or receiving signals
# num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
# weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
# gg <- list()
# for (i in 1:length(object.list)) {
#   gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax) + GG_KM_THEME
# }
# p <- cowplot::plot_grid(plotlist = gg, align = "hv", axis = "b")
# ggsave(paste0(figures_dir, plot_prefix, "_kos_re_compare_cell_pops.png"), plot = p, width = 10, height = 10)
# 
# 
# Identify the signaling changes of specific cell populations
# gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Mac") #, signaling.exclude = "MIF")
# gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "T") #, signaling.exclude = c("MIF"))
# 
# p <- patchwork::wrap_plots(plots = list(gg1,gg2))
# ggsave(paste0(figures_dir, plot_prefix, "_kos_re_compare_Mac_T.png"), plot = p)
# 
# # Identify signaling groups based on their functional similarity
# cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
# #> Compute signaling network similarity for datasets 1 2
# cellchat <- netEmbedding(cellchat, type = "functional")
# #> Manifold learning of the signaling networks for datasets 1 2
# cellchat <- netClustering(cellchat, type = "functional")
# #> Classification learning of the signaling networks for datasets 1 2
# # Visualization in 2D-space
# netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
# #> 2D visualization of signaling networks from datasets 1 2

# Identify altered signaling with distinct interaction strength
# (A) Compare the overall information flow of each signaling pathway or ligand-receptor pair
# CellChat can identify the conserved and context-specific signaling pathways by 
# simply comparing the information flow for each signaling pathway, which is defined by 
# the sum of communication probability among all pairs of cell groups in the inferred network (i.e., the total weights in the network).
signaling_contrib <- rankNet(cellchat, mode = "comparison", 
                             measure = "weight", 
                             sources.use = source_groups, 
                             targets.use = target_groups, 
                             return.data = T)$signaling.contribution %>% 
  pivot_wider(names_from = "group", values_from = c("contribution", "contribution.scaled", "contribution.relative.1")) %>% 
  mutate(scaled_diff = contribution.scaled_re - contribution.scaled_kos) %>% 
  arrange(desc(scaled_diff))

rankNet(cellchat, mode = "comparison", 
        measure = "weight", 
        signaling = signaling_contrib$name[1:10],
        sources.use = source_groups, 
        targets.use = target_groups, 
        stacked = F) + GG_KM_THEME

gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T) + GG_KM_THEME
gg2 <- rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F) + GG_KM_THEME
gg1 + gg2
ggsave(paste0(figures_dir, plot_prefix, "_inf_flow.png") , plot = gg1 + gg2, width = 14, height = 14)

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



# Identify dysfunctional signaling by comparing the communication probabities AND 
# the ligand-receptor differential gene expression and take the intersect of the two

# The above method for identifying the upgulated and down-regulated signaling is perfomed by comparing the communication probability between two datasets for each L-R pair and each pair of cell groups. 
# Alternatively, we can identify the upgulated and down-regulated signaling ligand-receptor pairs based on the differential expression analysis (DEA). Specifically, we perform differential expression analysis between two biological conditions (i.e., NL and LS) for each cell group, and then obtain the upgulated and down-regulated signaling based on the fold change of ligands in the sender cells and receptors in the receiver cells.
# Of note, users may observe the same LR pairs appearing in both the up-regulated and down-regulated results due to the fact that DEA between conditions is performed for each cell group. To perform DEA between conditions by ignoring cell group information, users can set group.DE.combined = TRUE in identifyOverExpressedGenes for CellChat v2.1.1.
# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pval_thresh <- 0.05
logFC_thresh <- 0.1
cellchat <- identifyOverExpressedGenes(cellchat, 
                                       group.dataset = "datasets", 
                                       pos.dataset = "re",
                                       only.pos = FALSE, 
                                       thresh.pc = 0.1, 
                                       thresh.fc = logFC_thresh, 
                                       thresh.p = pval_thresh, 
                                       group.DE.combined = FALSE)

# Map the results of differential expression analysis onto the inferred cell-cell 
# communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = "features", variable.all = TRUE)


combined_data <- c()
for (target in setdiff(levels(cellchat.kos@idents), "neuron")) {
  
  # RE L-R pairs with significant communication probabilities
  re.pairs.comprob <- netVisual_bubble(cellchat, 
                                       sources.use = source_groups, 
                                       targets.use = target_groups,  
                                       comparison = c(1, 2), 
                                       max.dataset = 2, 
                                       title.name = "Increased signaling in RE", 
                                       angle.x = 45, 
                                       thresh = pval_thresh, 
                                       sort.by.source = T, 
                                       remove.isolate = T, 
                                       return.data = T)$communication %>% 
    filter("re" == dataset) %>% 
    mutate(prob.scaled = scale(prob))
  
  # Extract the ligand-receptor pairs with upregulated ligands in RE
  re.pairs.dge <- subsetCommunication(cellchat, 
                                      net = net, 
                                      datasets = "re", 
                                      sources.use = source_groups, 
                                      targets.use = target_groups,   
                                      ligand.logFC = logFC_thresh)
  
  re_lr_intersect <- intersect(re.pairs.comprob$interaction_name, re.pairs.dge$interaction_name)
  
  re.pairs.comprob_and_dge <- re.pairs.comprob %>% 
    filter(interaction_name %in% re_lr_intersect)
    
  
  # KOS L-R pairs with significant communication probabilities
  kos.pairs.comprob <- netVisual_bubble(cellchat, 
                                        sources.use = source_groups, 
                                        targets.use = target_groups,  
                                        comparison = c(1, 2), 
                                        max.dataset = 1, 
                                        title.name = "Increased signaling in KOS", 
                                        angle.x = 45, 
                                        thresh = pval_thresh, 
                                        sort.by.source = T, 
                                        remove.isolate = T, 
                                        return.data = T)$communication %>% 
    filter("kos" == dataset) %>% 
    mutate(prob.scaled = scale(prob))
  
  # Extract the ligand-receptor pairs with upregulated ligands and upregulated receptors in KOS, i.e.,downregulated in RE
  kos.pairs.dge <- subsetCommunication(cellchat, 
                                       net = net, 
                                       datasets = "kos", 
                                       ligand.logFC = -logFC_thresh)
  
  kos_lr_intersect <- intersect(kos.pairs.comprob$interaction_name, kos.pairs.dge$interaction_name)
  
  kos.pairs.comprob_and_dge <- kos.pairs.comprob %>% 
    filter(interaction_name %in% kos_lr_intersect)
  
  combined_data <- bind_rows(combined_data, re.pairs.comprob_and_dge, kos.pairs.comprob_and_dge)
}


combined_data <- 
  combined_data %>% 
  mutate(cc = paste0(str_to_title(source), " -> ", str_to_title(target)), 
         lr_pair = interaction_name_2) 

combined_data %>% 
  ggplot(., aes(x = dataset, y = lr_pair, color = prob.scaled)) +
  geom_point(pch = 18, size = 7) +
  facet_wrap(~ cc, scales = "free", nrow = 1) + 
  theme_linedraw() + 
  guides(color = guide_colourbar(title = "Commun\nProb", ticks.colour = NA, frame.colour = "black")) + 
  scale_colour_gradientn(colors = rev(colorRampPalette(RColorBrewer::brewer.pal(n = 10, name = "Spectral"))(99)), 
                         na.value = "white", 
                         #limits=c(quantile(df$prob, 0,na.rm= T), quantile(df$prob, 1,na.rm= T)),
                         breaks = c(quantile(combined_data$prob.scaled, 0,na.rm= T), quantile(combined_data$prob.scaled, 1,na.rm= T)), 
                         labels = c("min","max")) +
  #guides(color = guide_colourbar(barwidth = 0.5, title = "Commun. Prob."))
  theme(legend.title = element_text(face = "bold", size = 10),
        #axis.text.x = element_tex(), 
        axis.text = element_text(face = "bold", size = 14),
        axis.title = element_blank(), 
        #aspect.ratio = 0.4, 
        panel.grid.major = element_blank()) + 
  geom_vline(xintercept=seq(1.5, length(unique(combined_data$cc)) -0.5, 1),lwd=0.1,colour="grey90") + 
  geom_hline(yintercept=seq(1.5, length(unique(combined_data$lr_pair)) -0.5, 1),lwd=0.1,colour="grey90")

## For TG: width = 14, height = 11
## For Cornea: width = , height = 
ggsave(paste0(figures_dir, plot_prefix, "_lr_pairs_up_down.png"), width = 14, height = 11)


##### PATHWAYS HEATMAP #####

## Set Source and Target Groups of Interest
if (grepl("tg", dataset, ignore.case = TRUE)) {
  
  # TG: interested in Neuron as source
  source_groups <- c("neuron")
  target_groups <- setdiff(levels(cellchat.kos@idents), c("neuron"))
  
} else if(grepl("cornea", dataset, ignore.case = TRUE)) {
  
  # Cornea: interested in Neuron as receiver
  source_groups <- setdiff(levels(cellchat.kos@idents), c("neuron"))
  target_groups <- c("neuron")
}

notable_pathways <- unique(combined_data$pathway_name)
df <- tibble()
conditions <- names(cellchat@netP)

for (condition in conditions) { 
  for (pathway in cellchat@netP[[condition]]$pathways) { 
    
    rows <- cellchat@netP[[condition]]$prob[, , pathway] %>% 
      as_tibble(rownames = "source") %>% 
      filter(source %in% source_groups) %>% 
      select(source, all_of(target_groups)) %>% 
      pivot_longer(all_of(target_groups), names_to = "target", values_to = "pathway_commprob") %>% 
      mutate(pathway = pathway, condition = condition)
    
    df <- bind_rows(df, rows)
  }
}

pathway_strength <- df %>% 
  group_by(source, target, pathway) %>% 
  filter(sum(pathway_commprob) > 0) %>%    ### Remove pathways that have a zero communication probability
  ungroup() %>% 
  group_by(source, target, condition) %>% 
  #mutate(pathway_commprob.scaled = scale(pathway_commprob)) %>% 
  mutate(pathway_commprob.scaled = pathway_commprob) %>% 
  ungroup()

pathway_strength %>% 
  mutate(cc = paste0(source, " -> ", target)) %>% 
  filter(pathway %in% unique(combined_data$pathway_name)) %>% 
  ggplot(., aes(x = condition, y = pathway, color = pathway_commprob.scaled)) +
  geom_point(pch = 18, size = 7) +
  facet_wrap(~ cc, scales = "free", nrow = 2) + 
  theme_linedraw() + 
  guides(color = guide_colourbar(title = "Commun\nProb", ticks.colour = NA, frame.colour = "black")) + 
  scale_colour_gradientn(colors = rev(colorRampPalette(RColorBrewer::brewer.pal(n = 10, name = "Spectral"))(99)), 
                         na.value = "white", 
                         #limits=c(quantile(df$prob, 0,na.rm= T), quantile(df$prob, 1,na.rm= T)),
                         breaks = c(quantile(combined_data$prob.scaled, 0,na.rm= T), quantile(combined_data$prob.scaled, 1,na.rm= T)), 
                         labels = c("min","max")) +
  #guides(color = guide_colourbar(barwidth = 0.5, title = "Commun. Prob."))
  theme(legend.title = element_text(face = "bold", size = 10),
        #axis.text.x = element_tex(), 
        axis.text = element_text(face = "bold", size = 14),
        axis.title = element_blank(), 
        #aspect.ratio = 0.4, 
        panel.grid.major = element_blank())

ggsave(paste0(figures_dir, plot_prefix, "_pathways_up_down.png"), width = 6, height = 9)



# Users can also find all the significant outgoing/incoming/both signaling according to the customized features and cell groups of interest
#df <- findEnrichedSignaling(object.list[[2]], features = c("CCL19", "CXCL12"), idents = c("Inflam. FIB", "COL11A1+ FIB"), pattern ="outgoing")




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

## Not sure what this one is
# weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
# par(mfrow = c(1,2), xpd=TRUE)
# for (i in 1:length(object.list)) {
#   netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
# }

#######################################################################################################################################

#######################################################################################################################################

# Extract the inferred cellular communication network as a data frame
# We provide a function subsetCommunication to easily access the inferred cell-cell communications of interest. For example,
# df.net <- subsetCommunication(cellchat) returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. Set slot.name = "netP" to access the the inferred communications at the level of signaling pathways
# df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5)) gives the inferred cell-cell communications sending from cell groups 1 and 2 to cell groups 4 and 5.
# df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb")) gives the inferred cell-cell communications mediated by signaling WNT and TGFb.




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