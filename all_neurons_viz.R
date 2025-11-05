library(tidyverse)
library(ggpubr)
library(ggprism)
library(Seurat)

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

##########################################

# Neurons - all conditions ----
all_neurons <- SeuratDisk::LoadH5Seurat(paste0(data_dir, "seurat/neurons_all_conditions.with_viral.h5Seurat"))

# control_neurons <- subset(all_neurons, subset = condition == "control") %>% SeuratDisk::SaveH5Seurat(paste0(data_dir, "seurat/control_neurons.h5Seurat"), overwrite = TRUE)

## Neuron Cluster Mean Normalized Expression Heatmaps ----
neuron_markers <- c("Trpv1", "Calca", "Tac1", "Gfra2", "Piezo2", "Trpm8", "P2rx3", "Mrgprd")
neuroimmune_markers <- c("Ccl2", "Cd44", "Tnfrsf1a", "Tnfrsf11a", "Tnfrsf21", "Vegfa", "Nrp1", "Nrp2", "Ly86", "Il10rb", 
                         "App", "Ccl21a", "Cd55", "Lgals9", "Mif", "Ptn")
renthal_neuron_atlas_markers <- c("Gfra2", "Pou4f2", "Tac1", "Gal", "Trpm8", "Cd55", "Scn11a", "Fxyd7", "Ngfr", "Nefh", 
                                  "Hapln4", "Cbln2", "Kcnab1", "Sst", "Il31ra", "Apoe", "Fabp7", "Mpz", "Gldn", "Scn7a", 
                                  "Dcn", "Pdgfra", "Mgp", "Alpl", "Cd74", "Igfbp7", "Tinagl1")

# UMAP by Condition
DimPlot(all_neurons, 
        group.by = "condition", 
        reduction = "umap", 
        pt.size = 1) + 
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

ggsave(paste0(figures_dir, "neurons_umap_by_condition.pdf"), width = 5, height = 5, bg = "white")

DimPlot(all_neurons, 
        group.by = "seurat_clusters",
        split.by = "condition",
        reduction = "umap", 
        pt.size = 1) + 
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

ggsave(paste0(figures_dir, "neurons_umap_by_cluster.pdf"), width = 20, height = 5, bg = "white")

DotPlot(all_neurons, 
        assay = "SCT", 
        features = renthal_neuron_atlas_markers, 
        #cols = viridis::viridis(n = 4),
        group.by = "seurat_clusters",
        #split.by = "condition",
        dot.scale = 6) + 
  #scale_fill_viridis_d() + 
  ggtitle("Neurons - All Conditions")

ggsave(paste0(figures_dir, "renthal_neuron_atlas_markers_dotplot.pdf"), width = 20, height = 10, bg = "white")

neuron_markers_heatmap <- 
  c("Trpm8", "Gfra2", "Mrgprd", "Trpv1", "Trpa1", "Tac1", "Calca", "Gfra3", "Ntrk1", "Ntrk2", "Ntrk3", "Ngfr", "Piezo2", "Nefh", "Ret")

immune_markers_heatmap <- 
  c("Tnfrsf1a", "Cd44", "Il13ra1", "Il10rb", "Vegfa", "Tnfrsf11a", "Tnfrsf21", "Hk1", "Ly86", "Npr2", "Il1rl1", "Il6ra", "Cd55", "Nrp1")

VlnPlot(all_neurons, 
        assay = "SCT",
        features = c(neuron_markers, neuroimmune_markers), 
        group.by = "seurat_clusters",
        split.by = "condition",
        stack = TRUE, 
        flip = TRUE) + 
  scale_fill_viridis_d() + 
  ggtitle("Neurons (All Conditions)")

ggsave(paste0(figures_dir, "neuron_markers_vln.png"), width = 10, height = 10, bg = "white")

VlnPlot(all_neurons, 
        assay = "SCT",
        features = c("Ccl2", neuron_markers_heatmap, immune_markers_heatmap), 
        group.by = "seurat_clusters",
        split.by = "condition",
        stack = TRUE, 
        flip = TRUE) + 
  scale_fill_viridis_d() + 
  ggtitle("Neurons (All Conditions)")

ggsave(paste0(figures_dir, "ccl2_neuron_and_immune_markers_vln.png"), width = 14, height = 14, bg = "white")

# CELL DEATH ----

neuronal_death_genes <- readr::read_csv("Mouse_Neuronal_Death_Injury_Gene_Panel.csv")
condition_avgs <- 
  AverageExpression(all_neurons, 
                    assays = "SCT", 
                    layer = "scale.data", 
                    group.by = "condition",
                    features = neuronal_death_genes$GeneSymbol)[[1]] %>%
  as.matrix() %>%
  as_tibble(rownames = "gene") %>%
  print(n = Inf)

condition_avgs

top_neuronal_death_genes <- c("Bax", "Bak1", "Bbc3", "Bcl2l11", "Bid", "Trp53", "Casp2", "Cycs", "Mlkl", "Zbp1", "Casp1", "Gsdmd", 
                              "Ddit3", "Oas1a", "Ifit1", "Ifit3", "Isg15", "Rsad2", "Sod2", "Atf4", "Ccl2", "Ccl5", 
                              "Cxcl10", "Stat1", "Irf7", "Gfap")

# Brief explanations
neuronal_death_genes %>% 
  filter(GeneSymbol %in% top_neuronal_death_genes) %>%
  print(n = Inf)

VlnPlot(all_neurons, 
        assay = "SCT",
        features = top_neuronal_death_genes, 
        group.by = "seurat_clusters",
        split.by = "condition",
        stack = TRUE, 
        flip = TRUE) + 
  scale_fill_viridis_d() + 
  ggtitle("Neurons (All Conditions)")

ggsave(paste0(figures_dir, "neuron_markers_vln.png"), width = 10, height = 10, bg = "white")

neuron_averages <- AverageExpression(all_neurons, 
                                     assays = "SCT", layer = "data",
                                     features = top_neuronal_death_genes,
                                     return.seurat = TRUE, 
                                     group.by = "condition")

DoHeatmap(neuron_averages, 
          features = top_neuronal_death_genes, 
          angle = 0, hjust = 0.5, vjust = 0,
          group.colors = rep("white", 4),
          raster = FALSE,
          draw.lines = FALSE) + 
  theme(axis.text = element_text(face = "bold", size = 18)) + 
  guides(color = "none") + 
  ggtitle("Back-labeled Corneal Sensory Neurons") + 
  guides(fill = guide_colourbar(title = "Mean Scaled Expression", 
                                title.theme = element_text(face = "bold", size = 14),
                                barheight = 1, 
                                barwidth = 20)) + 
  scale_y_discrete(position = "right") + 
  theme(legend.position = "bottom", 
        legend.title.position = "top",
        plot.title = element_text(face = "bold", hjust = 0, size = 14)) + 
  viridis::scale_fill_viridis()

ggsave(paste0(figures_dir, "neuron_death_genes_heatmap.png"), width = 5, height = 11, bg = "white")


# DotPlot(all_neurons, 
#         assay = "SCT",
#         features = c("Ccl2", neuron_markers_heatmap, immune_markers_heatmap), 
#         cols = viridis::viridis(n = 4),
#         group.by = "seurat_clusters",
#         split.by = "condition",
#         dot.scale = 6) + 
#   scale_fill_viridis_d() + 
#   ggtitle("Neurons - All Conditions")

# CONTROL NEURONS UMAP AND HEATMAP PLOTS ----


control_neurons <- all_neurons %>% subset(subset = condition == "control")

## DEG List
Idents(control_neurons) <- control_neurons$seurat_clusters
control_neurons %>% FindAllMarkers(assay = "SCT", ) %>% 
  select(cluster, gene, p_val, p_val_adj, everything()) %>% 
  group_by(cluster) %>% 
  arrange(p_val_adj, .by_group = TRUE) %>% 
  write_csv("table3_neuron_cluster_marker_genes.csv")


## Leiden clusters UMAP
neuron_umap <- 
  DimPlot(control_neurons, 
          group.by = "seurat_clusters", 
          reduction = "umap", 
          label = TRUE, 
          pt.size = 1.4,
          label.box = TRUE, 
          label.size = 5, repel = TRUE) + 
  NoLegend() + 
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))


#ggsave(paste0(figures_dir, "neuron_umap.svg"), width = 4, height = 4)


## Neuron Feature Plots
#FeaturePlot(control_neurons, features = neuron_markers, order = TRUE, pt.size = 1.4)

#library(patchwork)
feature_plots <- list()
feature_plots[["umap"]] <- neuron_umap + ggtitle(" ")  # To align with the feature plots below
for (i in neuron_markers) {
  feature_plots[[i]] <- 
    FeaturePlot(control_neurons, 
                features = i,
                #alphcontrol_neuronsa = 0.7,
                #cols = c("darkgreen"), 
                order = TRUE, 
                #max.cutoff = 1.0, 
                pt.size = 1.4) + NoLegend() + 
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"), 
          plot.title = element_text(face = "bold.italic", size = 14))
}

combined_feature_plots <- cowplot::plot_grid(plotlist = feature_plots, ncol = 3, axis = "bltr")
ggsave(paste0(figures_dir, "neuron_markers.svg"), plot = combined_feature_plots, width = 9, height = 9)


feature_plots <- list()
for (i in neuroimmune_markers) {
  feature_plots[[i]] <- 
    FeaturePlot(control_neurons, 
                features = i,
                #alphcontrol_neuronsa = 0.7,
                #cols = c("darkgreen"), 
                order = TRUE, 
                #max.cutoff = 1.0, 
                pt.size = 1.4) + NoLegend() + 
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"), 
          plot.title = element_text(face = "bold.italic", size = 14))
}

combined_feature_plots <- cowplot::plot_grid(plotlist = feature_plots, ncol = 4, axis = "bltr")
ggsave(paste0(figures_dir, "neuroimmune_markers.svg"), plot = combined_feature_plots, width = 12, height = 12)


## Neuron Markers violin plot
control_neurons@assays$SCT$data[c(neuron_markers, "Gfra3"), ] %>% 
  t() %>% 
  as_tibble(rownames = "barcode") %>% 
  pivot_longer(!barcode, names_to = "gene", values_to = "value") %>% 
  ggplot(aes(x = factor(gene, levels = c("Calca", "Trpv1", "Tac1", "Gfra3", "P2rx3", "Gfra2", "Mrgprd", "Trpm8", "Piezo2")), 
             y = value)) +
  geom_boxplot(size = 0.5) +
  geom_point(position = position_jitter(width = .1, seed = 0), size = 1, alpha = .5) + 
  ggprism::theme_prism() + 
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, vjust = 0.7, face = "bold.italic")) + 
  labs(y = "Normalized Expression")

ggsave(paste0(figures_dir, "neuron_markers_boxplot.svg"), width = 5, height = 5)


## Neuron Cluster Mean Normalized Expression Heatmaps ## 
neuron_markers_heatmap <- 
  c("Trpm8", "Gfra2", "Mrgprd", "Trpv1", "Trpa1", "Tac1", "Calca", "Gfra3", "Ntrk1", "Ntrk2", "Ntrk3", "Ngfr", "Piezo2", "Nefh", "Ret")

immune_markers_heatmap <- 
  c("Tnfrsf1a", "Cd44", "Il13ra1", "Il10rb", "Vegfa", "Tnfrsf11a", "Tnfrsf21", "Hk1", "Ly86", "Npr2", "Ccl2", "Il1rl1", "Il6ra", "Cd55", "Nrp1")


## Stacked plot - proportion of cells expressing each feature and color = intensity (i.e., mean residual of SCT)
scaled_neuron_counts <- control_neurons@assays$SCT$scale.data[neuron_markers_heatmap, ] %>% t() %>% 
  as_tibble(rownames = "barcode") %>% 
  inner_join(control_neurons@meta.data %>% select(seurat_clusters) %>% rownames_to_column("barcode"), 
             by = join_by("barcode" == "barcode")) %>% 
  group_by(seurat_clusters) %>% 
  summarise(across(!barcode, mean)) %>% 
  ungroup() %>% 
  column_to_rownames("seurat_clusters") %>% 
  as.matrix() %>% t()

heatmap <- 
  ComplexHeatmap::Heatmap(scaled_neuron_counts, 
                          rect_gp = gpar(col = "black", lwd = 2),
                          border = TRUE,
                          border_gp = gpar(col = "black", lwd = 2),
                          col = rev(RColorBrewer::brewer.pal(11,"RdBu")), 
                          cluster_rows = FALSE, 
                          cluster_columns = FALSE,
                          row_names_side = "left", 
                          row_names_gp = gpar(fontsize = 14, fontface = "bold.italic"),
                          column_names_gp = gpar(fontsize = 14, fontface = "bold"), 
                          column_names_side = "top",
                          column_names_rot = 0, 
                          column_names_centered = TRUE,
                          width = ncol(scaled_neuron_counts) * unit(7, "mm"), 
                          height = nrow(scaled_neuron_counts) * unit(7, "mm"), 
                          heatmap_legend_param = list(title = "Mean\nScaled\nExpr",
                                                      legend_height = unit(35, "mm"),
                                                      grid_width = unit(5, "mm"),
                                                      at = c(-5, 0, 5),
                                                      border = "black", just = "middle",
                                                      labels = c("-5", "0", "+5"))
  )

svg(paste0(figures_dir, "heatmap_nm.svg"), width = 5, height = 5)
heatmap
dev.off()


## Stacked plot - proportion of cells expressing each feature and color = intensity (i.e., mean residual of SCT)
scaled_neuron_counts <- control_neurons@assays$SCT$scale.data[immune_markers_heatmap, ] %>% t() %>% 
  as_tibble(rownames = "barcode") %>% 
  inner_join(control_neurons@meta.data %>% select(seurat_clusters) %>% rownames_to_column("barcode"), 
             by = join_by("barcode" == "barcode")) %>% 
  group_by(seurat_clusters) %>% 
  summarise(across(!barcode, mean)) %>% 
  ungroup() %>% 
  column_to_rownames("seurat_clusters") %>% 
  as.matrix() %>% t()

heatmap <- 
  ComplexHeatmap::Heatmap(scaled_neuron_counts, 
                          rect_gp = gpar(col = "black", lwd = 2),
                          border = TRUE,
                          border_gp = gpar(col = "black", lwd = 2),
                          col = rev(RColorBrewer::brewer.pal(11,"RdBu")), 
                          cluster_rows = FALSE, 
                          cluster_columns = FALSE,
                          row_names_side = "left", 
                          row_names_gp = gpar(fontsize = 14, fontface = "bold.italic"),
                          column_names_gp = gpar(fontsize = 14, fontface = "bold"), 
                          column_names_side = "top",
                          column_names_rot = 0, 
                          column_names_centered = TRUE,
                          width = ncol(scaled_neuron_counts) * unit(7, "mm"), 
                          height = nrow(scaled_neuron_counts) * unit(7, "mm"), 
                          heatmap_legend_param = list(title = "Mean\nScaled\nExpr",
                                                      legend_height = unit(35, "mm"),
                                                      grid_width = unit(5, "mm"),
                                                      at = c(-5, 0, 5),
                                                      border = "black", just = "middle",
                                                      labels = c("-5", "0", "+5"))
  )

svg(paste0(figures_dir, "heatmap_im.svg"), width = 5, height = 5)
heatmap
dev.off()





#######################################################################################################################################
## Generate TG Immune Umap with Markers

immune_markers <- c("Ptprc", "Pax5", "Cd19", "Ighm", 
                    "Cd3d", 
                    "Klrb1c", 
                    "Adgre1", "Cd83","Cd68")

combined_tg <- seurat_data <- SeuratDisk::LoadH5Seurat(paste0(data_dir, "seurat/combined_tg_annot.h5Seurat"))
# subset(subset = cell_L2 %in% c("Mye", "NK", "B", "T")) %>% 
# RunPCA(assay = "SCT", npcs = 50) %>% 
# RunUMAP(assay = "SCT", reduction = "pca", dims = 1:50)

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

ggsave(paste0(figures_dir, "tg_umap.svg"), width = 5, height = 5)

# combined_tg_props %>% 
#   ggplot(aes(cell_L2, prop)) + 
#   scale_fill_viridis_d() + 
#   geom_col(position=position_dodge()) + 
#   ggprism::theme_prism() + 
#   labs(x = "", y = "Proportion of TG immune cells") + 
#   theme(aspect.ratio = 1.7)

feature_plots <- list()
for (i in immune_markers) {
  feature_plots[[i]] <- 
    FeaturePlot(combined_tg, 
                features = i,
                #alpha = 0.7,
                #cols = c("darkgreen"), 
                order = TRUE, 
                raster = TRUE,
                max.cutoff = 1, 
                pt.size = 1.4) + NoLegend() + 
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"), 
          plot.title = element_text(face = "bold.italic"),
          axis.ticks = element_blank(), 
          axis.title = element_blank(), 
          axis.text = element_blank())
}

combined_tg_feature_plots <- cowplot::plot_grid(plotlist = feature_plots, ncol = 3, axis = "bltr")
ggsave(paste0(figures_dir, "tg_features.svg"), width = 7, height = 7)



## Generate Cornea Immune Umap with Markers
## Generate TG Immune Umap with Markers
sept_cornea <- seurat_data <- SeuratDisk::LoadH5Seurat(paste0(data_dir, "seurat/sept_cornea_annot.h5Seurat")) %>% 
  subset(subset = cell_L2 %in% c("Mye", "NK", "B", "T")) %>% 
  RunPCA(assay = "SCT", npcs = 50) %>% 
  RunUMAP(assay = "SCT", reduction = "pca", dims = 1:50)

# Immune Cell Proportions
sept_cornea_props <- sept_cornea@meta.data %>% 
  filter("control" == condition, cell_L2 %in% c("Mye", "NK", "B", "T")) %>% 
  group_by(condition, cell_L2) %>% 
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(prop = round(n / sum(n) * 100, digits = 1))

sept_cornea$umap_label <- factor(paste0(sept_cornea$cell_L2, "\n",
                                        sept_cornea_props$prop[as.numeric(sept_cornea$cell_L2)], 
                                        "%"))

cornea_umap <- DimPlot(sept_cornea, 
                       group.by = "umap_label", 
                       reduction = "umap", 
                       label = TRUE, label.box = TRUE, label.size = 5, repel = TRUE) + 
  NoLegend() + ggtitle("Cornea immune cell proportions")

ggsave(paste0(figures_dir, "cornea_umap.svg"), width = 5, height = 5)

feature_plots <- list()
for (i in immune_markers) {
  feature_plots[[i]] <- 
    FeaturePlot(sept_cornea, 
                features = i,
                #alpha = 0.7,
                #cols = c("darkgreen"), 
                order = TRUE, 
                raster = TRUE, 
                max.cutoff = 1, 
                pt.size = 1.4) + NoLegend() + 
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"), 
          plot.title = element_text(face = "bold.italic"),
          axis.ticks = element_blank(), 
          axis.title = element_blank(), 
          axis.text = element_blank())
}

combined_cornea_feature_plots <- cowplot::plot_grid(plotlist = feature_plots, ncol = 3, axis = "bltr")
ggsave(paste0(figures_dir, "cornea_features.svg"), width = 7, height = 7)


cowplot::plot_grid(plotlist = list(tg_umap, combined_tg_feature_plots, cornea_umap, combined_cornea_feature_plots), 
                   ncol = 2, axis = "bltr")

ggsave(paste0(figures_dir, "all_umap_features.svg"), width = 10, height = 10)


# FGF13 and related genes ----
#   |    | Gene                      | Why it’s closely related to FGF13 (and corneal pain)                                                                                                                                                             |
#   | -- | ------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
#   | 1  | **SCN9A (Nav1.7)**        | Principal sodium channel that FGF13 binds and tunes; its activity drives nociceptor excitability, and inhibiting Nav1.7 is analgesic. ([PMC][1], [MedlinePlus][2])                                               |
#   | 2  | **FGF14**                 | Fellow intracellular FGF (FHF) that cooperates with FGF13 to position Nav channels at the axon initial segment and modulate inactivation kinetics, influencing pain signalling. ([PNAS][3], [BioMed Central][4]) |
#   | 3  | **FGF12**                 | Another FHF; loss or mutation alters Nav1.2/1.6 inactivation, suggesting functional overlap or compensation for FGF13 in sensory neurons. ([BioMed Central][4], [PMC][5])                                        |
#   | 4  | **FGF11**                 | Fourth member of the intracrine FGF subfamily; shares structural domains with FGF13 and is predicted to modulate Nav gating in neurons. ([MDPI][6])                                                              |
#   | 5  | **ANK3 (Ankyrin-G)**      | Scaffold anchoring Nav-FGF complexes; required for proper clustering of sodium channels and therefore for FGF13-dependent modulation of excitability. ([PMC][7], [MDPI][8])                                      |
#   | 6  | **SPTBN4 (βIV-spectrin)** | Works with ankyrin-G to stabilize Nav channel clusters; disruption disperses Nav channels, indirectly impairing FGF13 function. ([PMC][9])                                                                       |
#   | 7  | **SCN8A (Nav1.6)**        | Sodium channel whose current amplitude/inactivation are enhanced by FGF13; contributes to action-potential propagation in sensory axons. ([Deep Blue][10], [BioMed Central][4])                                  |
#   | 8  | **SCN5A (Nav1.5)**        | Cardiac Nav isoform regulated by FGF13; illustrates a conserved FGF13–Nav interaction mechanism that can inform studies of other Navs in sensory neurons. ([PMC][11])                                            |
#   | 9  | **SCN2A (Nav1.2)**        | Nav1.2 gating is shifted by FHF proteins and likely sensitive to FGF13; expressed in some sensory neurons, offering another channel whose excitability can be shaped by FHFs. ([BioMed Central][4])              |
#   | 10 | **TRPV1**                 | Heat-activated TRP channel upstream of Nav1.7; TRPV1-mediated heat stimuli engage the FGF13–Nav1.7 pathway during nociception, linking this cation channel to corneal pain modulation. ([PMC][1])                |
#   
#   [1]: https://pmc.ncbi.nlm.nih.gov/articles/PMC12259270/?utm_source=chatgpt.com "Sensory neuron–expressed FGF13 controls nociceptive signaling in ..."
# [2]: https://medlineplus.gov/genetics/gene/scn9a/?utm_source=chatgpt.com "SCN9A gene: MedlinePlus Genetics"
# [3]: https://www.pnas.org/doi/10.1073/pnas.1521194113?utm_source=chatgpt.com "Polarized localization of voltage-gated Na+ channels is ... - PNAS"
# [4]: https://biosignaling.biomedcentral.com/articles/10.1186/s12964-023-01284-0?utm_source=chatgpt.com "Fibroblast growth factor signaling in axons: from development to ..."
# [5]: https://pmc.ncbi.nlm.nih.gov/articles/PMC9429545/?utm_source=chatgpt.com "Modulating effects of FGF12 variants on NaV1.2 and NaV1.6 being ..."
# [6]: https://www.mdpi.com/2813-3137/3/1/5?utm_source=chatgpt.com "Pathology and Therapeutic Significance of Fibroblast Growth Factors"
# [7]: https://pmc.ncbi.nlm.nih.gov/articles/PMC2133082/?utm_source=chatgpt.com "AnkyrinG Is Required for Clustering of Voltage-gated Na Channels ..."
# [8]: https://www.mdpi.com/2218-273X/15/6/901?utm_source=chatgpt.com "Ankyrin-G and Its Binding Partners in Neurons - MDPI"
# [9]: https://pmc.ncbi.nlm.nih.gov/articles/PMC2199236/ "
#             βIV-spectrin regulates sodium channel clustering through ankyrin-G at axon initial segments and nodes of Ranvier - PMC
#         "
# [10]: https://deepblue.lib.umich.edu/bitstream/2027.42/100063/1/janelleo_1.pdf?utm_source=chatgpt.com "[PDF] Regulation and Mutation of Voltage-Gated Sodium Channel SCN8A ..."
# [11]: https://pmc.ncbi.nlm.nih.gov/articles/PMC3383600/?utm_source=chatgpt.com "Fibroblast growth factor homologous factor 13 regulates Na+ ..."

fgf13_genes <- c("Fgf13", "Scn9a", "Fgf14", "Fgf12", "Fgf11", "Ank3", "Sptbn4", "Scn8a", "Scn5a", "Scn2a", "Trpv1")
neuron_markers <- c("Trpv1", "Calca", "Tac1", "Gfra2", "Piezo2", "Trpm8", "P2rx3", "Mrgprd")

VlnPlot(all_neurons, 
        assay = "SCT",
        features = c(fgf13_genes), 
        #group.by = "seurat_clusters",
        group.by = "condition",
        #split.by = "condition",
        stack = TRUE, 
        flip = TRUE) + 
  scale_fill_viridis_d() + 
  ggtitle("Back-labeled Corneal Sensory Neurons")

ggsave(paste0(figures_dir, "neuron_fgf13_genes_vln.png"), width = 5, height = 5, bg = "white")

neuron_averages <- AverageExpression(all_neurons, 
                                     assay = "SCT", layer = "data",
                                     features = fgf13_genes,
                                     return.seurat = TRUE, 
                                     group.by = "condition")

DoHeatmap(neuron_averages, 
          features = fgf13_genes, 
          angle = 0, hjust = 0.5, vjust = 0,
          group.colors = rep("white", 4),
          raster = FALSE,
          draw.lines = FALSE) + 
  theme(axis.text = element_text(face = "bold", size = 18)) + 
  guides(color = "none") + 
  ggtitle("Back-labeled Corneal Sensory Neurons") + 
  guides(fill = guide_colourbar(title = "Mean Scaled Expression", 
                                title.theme = element_text(face = "bold", size = 14),
                                barheight = 1, 
                                barwidth = 20)) + 
  scale_y_discrete(position = "right") + 
  theme(legend.position = "bottom", 
        legend.title.position = "top",
        plot.title = element_text(face = "bold", hjust = 0, size = 14)) + 
  viridis::scale_fill_viridis()

ggsave(paste0(figures_dir, "neuron_fgf13_genes_heatmap.png"), width = 5, height = 5, bg = "white")

# Viral Genes ----

kos_p <- 
  VlnPlot(all_neurons %>% subset(subset = condition == "kos"), 
        assay = "SCT",
        features = c(neuron_markers, neuroimmune_markers), 
        group.by = "seurat_clusters",
        split.by = "viral_infection",
        stack = TRUE, 
        flip = TRUE) + 
  scale_fill_viridis_d() + 
  ggtitle("Neurons (KOS)")

re_p <-
  VlnPlot(all_neurons %>% subset(subset = condition == "re"), 
        assay = "SCT",
        features = c(neuron_markers, neuroimmune_markers), 
        group.by = "seurat_clusters",
        split.by = "viral_infection",
        stack = TRUE, 
        flip = TRUE) + 
  scale_fill_viridis_d() + 
  ggtitle("Neurons (RE)")

# Make a combined plot with kos vertically on top of re
cowplot::plot_grid(kos_p, re_p, ncol = 1, align = "v", axis = "lr")

ggsave(paste0(figures_dir, "infected_violin.png"), width = 10, height = 20, bg = "white")


# Percent Infected ----
key_markers <- c(neuron_markers_heatmap, immune_markers_heatmap)

# load the counts from the neurons with the viral transcripts (to get our original clustering)
all_neurons.no_viral <- SeuratDisk::LoadH5Seurat(paste0(data_dir, "seurat/neurons_all_conditions.h5Seurat"))

# Create a dataframe with barcodes and their original clusters from all_neurons.no_viral
original_clusters <- all_neurons.no_viral@meta.data %>%
  select(seurat_clusters) %>%
  rename(original_clusters = seurat_clusters) %>%
  rownames_to_column("barcode") %>% 
  mutate(barcode = paste0("neuron_", barcode)) %>% 
  column_to_rownames("barcode")

# Set a new metadata field "original_cluster" on the all_neurons object based on matching barcodes with original_clusters
all_neurons <- AddMetaData(
  object = all_neurons,
  metadata = original_clusters,
  col.name = 'original_clusters'
)

# Barplot of proportion of infected cells in each original cluster, faceted by condition
percent_infected_by_cluster <-
  all_neurons@meta.data %>% 
  filter(condition %in% c("kos", "re")) %>%
  select(original_clusters, condition, viral_infection) %>% 
  ggplot(aes(x = condition, fill = viral_infection)) + 
  geom_bar(position = "fill") + 
  facet_wrap(~ original_clusters, ncol = length(unique(all_neurons@meta.data$original_clusters))) + 
  scale_fill_prism(palette = "colorblind_safe") + 
  ggpubr::theme_pubr() + 
  labs(y = "Proportion of Cells\nin Cluster for Condition", fill = "Viral Infection") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.7, size = 8),
        axis.title.x = element_blank(),
        strip.text = element_text(face = "bold", size = 6),
        axis.title.y = element_text(face = "bold", size = 8),
        axis.text.y = element_text(size = 8),
        legend.title = element_text(face = "bold", size = 8),
        legend.text = element_text(size = 6))

ggsave(paste0(figures_dir, "percent_infected_by_cluster.png"), plot = percent_infected_by_cluster, width = 4, height = 2, bg = "white")

# Heatmaps
scaled_neuron_counts <- all_neurons@assays$SCT$scale.data[key_markers, ] %>% t() %>% 
    as_tibble(rownames = "barcode") %>% 
    inner_join(all_neurons@meta.data %>% select(original_clusters) %>% rownames_to_column("barcode"), 
               by = join_by("barcode" == "barcode")) %>% 
    group_by(original_clusters) %>% 
    summarise(across(!barcode, mean)) %>% 
    ungroup() %>% 
    column_to_rownames("original_clusters") %>% 
    as.matrix() %>% t()
  
  heatmap <- 
    ComplexHeatmap::Heatmap(scaled_neuron_counts, 
                            rect_gp = grid::gpar(col = "black", lwd = 2),
                            border = TRUE,
                            border_gp = grid::gpar(col = "black", lwd = 2),
                            col = rev(RColorBrewer::brewer.pal(11,"RdBu")), 
                            cluster_rows = FALSE, 
                            cluster_columns = FALSE,
                            row_names_side = "left", 
                            row_names_gp = grid::gpar(fontsize = 14, fontface = "bold.italic"),
                            column_names_gp = grid::gpar(fontsize = 14, fontface = "bold"), 
                            column_names_side = "top",
                            column_names_rot = 0, 
                            column_names_centered = TRUE,
                            width = ncol(scaled_neuron_counts) * unit(7, "mm"), 
                            height = nrow(scaled_neuron_counts) * unit(7, "mm"), 
                            heatmap_legend_param = list(title = "Mean\nScaled\nExpr",
                                                        legend_height = unit(35, "mm"),
                                                        grid_width = unit(5, "mm"),
                                                        at = c(-5, 0, 5),
                                                        border = "black", just = "middle",
                                                        labels = c("-5", "0", "+5"))
    )
  
png(paste0(figures_dir, "heatmap_markers.png"), width = 5, height = 10, units = "in", res = 300)
heatmap
dev.off()
  