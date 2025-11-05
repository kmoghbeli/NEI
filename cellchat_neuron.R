library(tidyverse)
library(ggprism)
library(Seurat)
library(SeuratDisk)
library(CellChat)
library(patchwork)
library(ComplexHeatmap)
library(network)
library(sna)
library(ggnetwork)
options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 1000*1024^2 * 16)
conflicted::conflict_prefer_all("dplyr", quiet = TRUE)

config_dir <- "./config/"
data_dir <- "../data_objects/"
figures_dir <- "../figures/"

# Parallell Processing
if(parallel::detectCores(logical=FALSE) > 3) {
  library(doParallel)
  
  num_cores <- parallel::detectCores(logical=FALSE)
  cl <- makePSOCKcluster(num_cores - 2)
  registerDoParallel(cl)
}

# Set default font to Arial
theme_update(text = element_text(family = 'Arial'))

#######################################################################################################################################
##### PLOTTING FUNCTIONS #####

## Celltype-to-Celltype Circle Plot
cc_circle_plot <- function(cellchat_obj, source_groups, target_groups) {
  
  # Perform Max scaling of matrix subset of source/target groups of interest
  cc_mtx <- 
    cellchat_obj@net$weight %>% 
    as.data.frame() %>% 
    rownames_to_column("source") %>% 
    filter(source %in% levels(cellchat_obj@idents)[source_groups]) %>% 
    column_to_rownames("source") %>% 
    select(levels(cellchat_obj@idents)[target_groups])
  
  cc_mtx %>% 
    mutate(across(everything(), ~ (.x - min(cc_mtx)) / (max(cc_mtx) - min(cc_mtx)))) %>% 
    rownames_to_column("source") %>% 
    pivot_longer(!source, names_to = "target", values_to = "strength") %>% 
    ggnetwork(layout = "circle", arrow.gap = 0.055) %>% 
    ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_edges(aes(color = strength), size = 2, 
               arrow = arrow(length = unit(14, "pt"), type = "closed")) +
    #geom_nodes(color = alpha("black", 1), shape = 15) +  #aes(color = family, size = importance)) + 
    geom_nodelabel(aes(label = vertex.names), fontface = "bold", color = "black", size = 5, label.padding = unit(0.2, "lines"),) + 
    #geom_nodelabel(aes(label = vertex.names), fontface = "bold", color = "black", size = 5, label.padding = unit(0.25, "lines"),) + 
    #geom_edgelabel(aes(label = round(strength, 2)), size = 6, fill = alpha("white", 1)) + 
    guides(color = guide_colourbar(title = "Relative Communication Strength", title.position = "top", vjust = 0.5, ticks.colour = NA, frame.colour = "black")) + 
    scale_colour_gradientn(colors = colorRampPalette(RColorBrewer::brewer.pal(n = 10, name = "Greys")[3:9])(10), 
                           #na.value = "white", 
                           #limits=c(quantile(df$prob, 0,na.rm= T), quantile(df$prob, 1,na.rm= T)),
                           breaks = c(0, 1), 
                           labels = c("Low", "High")
    ) +
    theme_blank() + 
    #coord_cartesian(xlim=c(-0.5, 1.5)) + 
    scale_x_continuous(limits = c(-0.25, 1.25)) + 
    scale_y_continuous(limits = c(-0.1, 1.1)) + 
    theme(legend.position = "bottom", 
          legend.title.align = 0.5, 
          plot.title = element_text(face = "bold", size = 12, hjust = 0.5, vjust = 0.5),
          aspect.ratio = 1.1)
}

## PATHWAYS HEATMAP ##
pathways_heatmap <- function(cellchat_obj, signaling_poi, source_groups, target_groups, title, show_legend = TRUE) {

  source <- levels(cellchat_obj@idents)[source_groups]
  target <- levels(cellchat_obj@idents)[target_groups]
  
  df <- data.frame()
  
  for (pathway in cellchat_obj@netP$pathways) {
    row <- cellchat_obj@netP$prob[, , pathway][source, target] %>% 
      t() %>% 
      as.data.frame() %>% 
      mutate(pathway = pathway)
    
    df <- bind_rows(df, row)
  }
  
  pathway_relative_strength <- df %>% column_to_rownames("pathway")
  
  colnames(pathway_relative_strength) <- paste0(source, " -> ", target) %>% str_replace("neuron", "Neuron")
  
  # Remove rows for pathways with all zeroes
  pathway_relative_strength <- pathway_relative_strength %>% filter_all(any_vars(. != 0))
  
  pathway_relative_strength <- pathway_relative_strength / max(pathway_relative_strength)
  
  pathway_relative_strength.scaled <- scale(pathway_relative_strength)
  
  ## Filter down to pathways with a scaled strength (z-score) > threshold i.e., 0.8
  # pathway_z_score_threshold <- 1.3  # z-score 1.3 ~ 90%
  # 
  # pathway_relative_strength.scaled.significant <-
  #   pathway_relative_strength.scaled  %>%
  #   as.data.frame() %>%
  #   filter(if_any(everything(), ~ . > 1.0)) %>%
  #   as.matrix()
  # 
  # signaling_poi <- rownames(pathway_relative_strength.scaled.significant)   ## will use this in next step to filter the L-R pairs to our Pathways of interest (POI)
  
  na_pathways <- setdiff(signaling_poi, rownames(pathway_relative_strength.scaled))
  na_vals <- data.frame(matrix(data = NA, ncol = ncol(pathway_relative_strength.scaled), nrow = length(na_pathways)))
  rownames(na_vals) <- na_pathways
  colnames(na_vals) <- colnames(pathway_relative_strength.scaled)
  
  mtx_min <- min(pathway_relative_strength.scaled)
  mtx_max <- max(pathway_relative_strength.scaled)
  
  pathway_relative_strength.scaled <- rbind(pathway_relative_strength.scaled, na_vals)
  
  heatmap_mtx <- pathway_relative_strength.scaled[intersect(signaling_poi, rownames(pathway_relative_strength.scaled)), ] %>% as.matrix()
  
  ComplexHeatmap::Heatmap(heatmap_mtx,
                          col = rev(RColorBrewer::brewer.pal(11,"Spectral")), 
                          column_title = title, 
                          column_title_side = "top", column_title_gp = gpar(fontsize = 12, fontface = "bold"),
                          cluster_rows = FALSE, 
                          cluster_columns = FALSE,
                          row_names_side = "left", 
                          row_title = "Pathways", row_title_gp = gpar(fontsize = 18, fontface = "bold"),
                          row_names_gp = gpar(fontsize = 14, fontface = "bold"),
                          column_names_gp = gpar(fontsize = 14, fontface = "bold"), 
                          column_names_rot = 45, 
                          width = ncol(heatmap_mtx) * unit(10, "mm"), 
                          height = nrow(heatmap_mtx) * unit(7, "mm"), 
                          show_heatmap_legend = show_legend,
                          heatmap_legend_param = list(title = "Commun\nProb",
                                                      at = c(mtx_min, mtx_max),
                                                      legend_height = unit(35, "mm"),
                                                      grid_width = unit(5, "mm"), 
                                                      border = "black", just = "middle",
                                                      labels = c("min", "max")))
}

#lr_bubble_plot(cellchat.tg, neuron_rec_signaling_poi, tg_immune, tg_neuron, "TG")

## LIGAND-RECEPTOR BUBBLE PLOTS
lr_bubble_plot <- function(cellchat_obj, signaling_poi, source_groups, target_groups, title, show_legend = TRUE) {
  
  lr_pairs <- 
    netVisual_bubble(cellchat_obj, 
                     sources.use = source_groups, 
                     targets.use = target_groups, 
                     color.heatmap = c("Spectral"), 
                     sort.by.source.priority = T,
                     thresh = 0.05, 
                     font.size = 14, 
                     dot.size.min = 7, 
                     return.data = TRUE)$communication

  sigLRs <- cellchat_obj@LR$LRsig %>% filter(pathway_name %in% signaling_poi) %>% pull(interaction_name_2)
  
  lr_pairs.scaled <- lr_pairs %>% 
    mutate(label = paste0(source, " -> ", target), 
           lr_pair = interaction_name_2) %>% 
    select(label, prob, lr_pair) %>% 
    pivot_wider(names_from = "label", values_from = "prob") %>% 
    mutate(across(!starts_with("lr_pair"), scale)) %>% 
    filter(lr_pair %in% sigLRs) %>%   # NOTE: Filter AFTER scaling
    filter(if_any(!lr_pair, ~ . >= 0.5)) %>%   # NOTE: Filter AFTER scaling to keep only L-R pairs with a scaled value above our threshold
    pivot_longer(!lr_pair, names_to = "cc", values_to = "comm_prob") %>% 
    mutate(cc = str_replace(cc, "neuron", "Neuron")) %>% 
    drop_na(comm_prob)

  lr_pairs.scaled %>% 
    ggplot(., aes(x = factor(cc, levels = paste0(levels(cellchat_obj@idents)[source_groups], " -> ", levels(cellchat_obj@idents)[target_groups])), 
                  y = lr_pair, color = comm_prob)) +
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
          #legend.key.size = unit(100, 'mm'),
          axis.text = element_text(face = "bold", size = 18),
          axis.text.x = element_text(angle = 45, hjust = 1), 
          axis.text.y = element_text(face = "bold.italic"), 
          axis.title = element_blank(), 
          #aspect.ratio = 0.4, 
          panel.grid.major = element_blank(), 
          plot.title = element_text(face = "bold", size = 32, hjust = 0.5, vjust = 0.5)
    ) + 
    ggtitle(title) + 
    geom_vline(xintercept=seq(1.5, length(unique(lr_pairs.scaled$lr_pair)) -0.5, 1),lwd=0.1,colour="grey90") + 
    geom_hline(yintercept=seq(1.5, length(unique(lr_pairs.scaled$cc)) -0.5, 1),lwd=0.1,colour="grey90")
}

#lr_bubble_plot(cellchat.tg, neuron_send_signaling_poi, tg_neuron, tg_immune, "TG")


#######################################################################################################################################

cellchat.tg <- readr::read_rds(paste0(data_dir, "cellchat/cellchat_tg_control_with_neurons_celltypeloc.rds"))
cellchat.cornea <- readr::read_rds(paste0(data_dir, "cellchat/cellchat_cornea_control_with_neurons_celltypeloc.rds"))

tg_neuron <- match(c("Neuron"), levels(cellchat.tg@idents))
tg_immune <- match(setdiff(levels(cellchat.tg@idents), c("Neuron")), 
                   levels(cellchat.tg@idents))

cornea_neuron <- match(c("Neuron"), levels(cellchat.cornea@idents))
cornea_immune <- match(setdiff(levels(cellchat.cornea@idents), c("Neuron")), 
                   levels(cellchat.tg@idents))


## Cell-Cell Plots
tg_forward_cc <- cc_circle_plot(cellchat.tg, tg_neuron, tg_immune) + ggtitle("TG") + theme(legend.position="none")
cornea_forward_cc <- cc_circle_plot(cellchat.cornea, cornea_neuron, cornea_immune) + ggtitle("Cornea") + theme(legend.position="none")

tg_forward_cc + cornea_forward_cc + patchwork::plot_layout(guides = "collect", ncol = 2) &
  theme(legend.position = "bottom", legend.direction = "horizontal")

ggsave(paste0(figures_dir, "circle_plot_forward.svg"), width = 5, height = 5)

tg_reverse_cc <- cc_circle_plot(cellchat.tg, tg_immune, tg_neuron) + ggtitle("TG")
cornea_reverse_cc <- cc_circle_plot(cellchat.cornea, cornea_immune, cornea_neuron) + ggtitle("Cornea")

tg_reverse_cc + cornea_reverse_cc + patchwork::plot_layout(guides = "collect", ncol = 2) &
  theme(legend.position = "bottom", legend.direction = "horizontal")

ggsave(paste0(figures_dir, "circle_plot_reverse.svg"), width = 5, height = 5)



## Pathways Heatmaps
neuron_send_signaling_poi <- c("LAMININ", "COLLAGEN", "GALECTIN", "PTN", "APP", "MIF", "Cholesterol", "ADGRL", "ADGRE", "CCL", "CSF", "CypA", "CD39") #  "PARs", 
neuron_rec_signaling_poi <- c("LAMININ", "APP", "CNTN", "CypA", "Prostaglandin", "FN1", "CCL", "THBS", "TGFb", "PSAP", "CXCL", "PARs")

tg_forward_heatmap <- pathways_heatmap(cellchat.tg, neuron_send_signaling_poi, tg_neuron, tg_immune, "TG")
cornea_forward_heatmap <- pathways_heatmap(cellchat.cornea, neuron_send_signaling_poi, cornea_neuron, cornea_immune, "Cornea", show_legend = FALSE)

svg(paste0(figures_dir, "pathway_heatmap_forward.svg"), width = 14, height = 7)
draw(tg_forward_heatmap + cornea_forward_heatmap, ht_gap = unit(1, "cm"))
dev.off()


tg_reverse_heatmap <- pathways_heatmap(cellchat.tg, neuron_rec_signaling_poi, tg_immune, tg_neuron, "TG")
cornea_reverse_heatmap <- pathways_heatmap(cellchat.cornea, neuron_rec_signaling_poi, cornea_immune, cornea_neuron, "Cornea", show_legend = FALSE)

svg(paste0(figures_dir, "pathway_heatmap_reverse.svg"), width = 14, height = 7)
draw(tg_reverse_heatmap + cornea_reverse_heatmap, ht_gap = unit(1, "cm"))
dev.off()


tg_forward_bubble <- lr_bubble_plot(cellchat.tg, neuron_send_signaling_poi, tg_neuron, tg_immune, "TG")
cornea_forward_bubble <- lr_bubble_plot(cellchat.cornea, neuron_send_signaling_poi, cornea_neuron, cornea_immune, "Cornea", show_legend = FALSE)

tg_forward_bubble + cornea_forward_bubble + patchwork::plot_layout(guides = "collect", ncol = 2) &
  theme(legend.position = "right", 
        legend.direction = "vertical", 
        legend.key.size = unit(2, 'cm'), 
        legend.title = element_text(size=25, face = "bold"),
        legend.text = element_text(size=20))

ggsave(paste0(figures_dir, "bubbles_forward.svg"), width = 12, height = 27)

tg_reverse_bubble <- lr_bubble_plot(cellchat.tg, neuron_rec_signaling_poi, tg_immune, tg_neuron, "TG")
cornea_reverse_bubble <- lr_bubble_plot(cellchat.cornea, neuron_rec_signaling_poi, cornea_immune, cornea_neuron, "Cornea")

tg_reverse_bubble + cornea_reverse_bubble + patchwork::plot_layout(guides = "collect", ncol = 2) &
  theme(legend.position = "right", 
        legend.direction = "vertical", 
        legend.key.size = unit(1.5, 'cm'), 
        legend.title = element_text(size=25, face = "bold"),
        legend.text = element_text(size=20))

ggsave(paste0(figures_dir, "bubbles_reverse.svg"), width = 12, height = 27)


#######################################################################################################################################
###### BUILD CELLCHAT OBJECT [ONLY NEED TO RUN THIS ONCE TO PREP THE DATASET] ######

location_subset <- "TG"
dataset <- paste0(str_to_lower(location_subset), "_control_with_neurons")
grouping <- "celltypeloc"

plot_prefix <- paste0(dataset, "_", grouping, "_", ifelse(neuron_as_sender, "neuron_send_", "neuron_rec_"))

## 1) Load the full annotated, normalized dataset  
seurat_data <- SeuratDisk::LoadH5Seurat(paste0(data_dir, "seurat/control_with_neurons.h5Seurat"))  

seurat_data <- subset(seurat_data, subset = location %in% c(str_to_lower(location_subset), "corneal_afferents") & cell_L2 %in% c("neuron", "Mye", "NK", "B", "T"))

# Re-level celltypeloc
seurat_data$celltypeloc <- factor(seurat_data$cell_L2, levels = c("neuron", "Mye", "NK", "B", "T"), labels = c("Neuron", "Mye", "NK", "B", "T"))

# CellChat likes to have this metadata
seurat_data$samples <- factor(paste0(seurat_data$location, "_", seurat_data$date, "_", seurat_data$condition))

cellchat <- createCellChat(object = seurat_data, group.by = "celltypeloc", assay = "SCT")

cellchat@DB <- CellChatDB.mouse # use CellChatDB.human if running on human data

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

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

