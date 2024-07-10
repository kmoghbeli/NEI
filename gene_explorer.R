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
combined_tg <- SeuratDisk::LoadH5Seurat(paste0(data_dir, "combined_tg.h5Seurat"))

#sept_cornea <- SeuratDisk::LoadH5Seurat(paste0(data_dir, "sept_cornea_sct.h5Seurat"))
sept_cornea <- SeuratDisk::LoadH5Seurat(paste0(data_dir, "sept_cornea.h5Seurat"))
july_cornea <- SeuratDisk::LoadH5Seurat(paste0(data_dir, "july_cornea_sct.h5Seurat"))


##########################################

# We had a list of some genes we would like to see. 
# Specifically, we are hypothesizing that these genes would be in the epithelial cells at the cornea. 
# Epithelial cells, if I remember correctly from the UMAPs would be those that are expressing K14, K15 or other keratins. 
# Within that populations, we are wondering what the expression of these genes are:

epi_genes_of_interest <- c("Epcam", "Ptprc", "Ngf", "Bdnf", "Gdnf", "Artn", "Nrtn", "Ntf3", "Ntf5")

# Additionally, we would like to know if immune cells within the TG are expressing factors like IL-10, CCR2, and VEGF.
tg_immune_genes <- c("Il10", "Ccr2", "Vegfa")

sept_cornea.epi <- subset(sept_cornea, subset = "epithelial" == cell_L1)
july_cornea.epi <- subset(july_cornea, subset = "epithelial" == cell_L1)

VlnPlot(sept_cornea.epi, assay = "SCT",
        features = epi_genes_of_interest,
        group.by = "seurat_clusters", split.by = "condition",
        scCustomize::DiscretePalette_scCustomize(4, palette = "glasbey", shuffle = FALSE),
        stack = TRUE, flip = TRUE) + 
  ggtitle("Sept Cornea Epithelial Cells")

ggsave(paste0(figures_dir, "sept_cornea_epi_neuro_genes.png"), height = 7, width = 7, bg = "white")

VlnPlot(sept_cornea, assay = "SCT",
        features = epi_genes_of_interest,
        group.by = "seurat_clusters", split.by = "condition",
        scCustomize::DiscretePalette_scCustomize(4, palette = "glasbey", shuffle = FALSE),
        stack = TRUE, flip = TRUE) + 
  ggtitle("Sept Cornea All Cells")

ggsave(paste0(figures_dir, "sept_cornea_all_neuro_genes.png"), height = 14, width = 14, bg = "white")


VlnPlot(july_cornea.epi, assay = "SCT",
        features = epi_genes_of_interest,
        group.by = "seurat_clusters", split.by = "condition",
        scCustomize::DiscretePalette_scCustomize(4, palette = "glasbey", shuffle = FALSE),
        stack = TRUE, flip = TRUE) + 
  ggtitle("July Cornea Epithelial Cells")

ggsave(paste0(figures_dir, "july_cornea_epi_neuro_genes.png"), height = 7, width = 7, bg = "white")

VlnPlot(combined_tg, assay = "SCT",
        features = tg_immune_genes,
        group.by = "cell_L2", split.by = "condition",
        scCustomize::DiscretePalette_scCustomize(4, palette = "glasbey", shuffle = FALSE),
        stack = TRUE, flip = TRUE) + 
  ggtitle("June + August TG Immune Cells")

ggsave(paste0(figures_dir, "combined_tg_immune_genes.png"), height = 7, width = 7, bg = "white")

##########################################

# brian_genes <- c("Ccr2", 
#                  "Ifngr1", "Ifngr2", "Ifnar1", "Ifnar2", 
#                  "Tnfaip1", "Tnfaip3", "Tnfrsf21", "Tnfrsf11a", "Tnfrsf1a", 
#                  "Il4ra", "Il6", "Il6st", "Il10rb", "Il13ra1"
# )
# nt_receptor_genes <- c("Nmur1", "Adrb2", 
#                        #"Calcr", 
#                        "Vipr2", "Chrm1", 
#                        #"Chrna1", "Chrna7", 
#                        "Chrnb1", "Fpr1", 
#                        #"Mrgpra1", 
#                        "Mrgprb1", "Grm5", "Grm7", "Grin1", "Grin2a", "Tacr1", 
#                        "Bmpr2", "Csf3r")

