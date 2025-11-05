library(tidyverse)
library(ggpubr)
library(ggprism)
library(Seurat)
conflicted::conflict_prefer_all("dplyr", quiet = TRUE)

config_dir <- "./config/"
data_dir <- "../data_objects/"

# Parallell Processing
if (parallel::detectCores(logical = FALSE) > 3) {
  library(doParallel)

  num_cores <- parallel::detectCores(logical = FALSE)
  cl <- makePSOCKcluster(num_cores - 2)
  registerDoParallel(cl)
}

#######################################################################################################################################

dataset <- "combined_tg_viral_annot"

# Combine SEACells annotations with Seurat object ----


# 1) Load the full annotated, normalized dataset

seurat_data <- SeuratDisk::LoadH5Seurat(paste0(data_dir, "seurat/", dataset, ".h5Seurat"))

# 2) Load the SEAcells assignments

seacells_assignments <-
  read_csv(paste0(data_dir, "seacells/", dataset, ".seacells.assignments.csv")) %>%
  column_to_rownames("index")

# 3) Add the SEAcells assignments to the Seurat object

seurat_data <- AddMetaData(seurat_data, seacells_assignments)

seurat_data$cell_L2 <- factor(seurat_data$cell_L2)
seurat_data$SEACell <- factor(seurat_data$SEACell)

# 4) Save the Seurat object

SeuratDisk::SaveH5Seurat(seurat_data, paste0(data_dir, "seurat/", dataset, ".with_seacells.h5Seurat"), overwrite = TRUE)

# 5) Delete the old annotated dataset

file.remove(paste0(data_dir, "seurat/", dataset, ".h5Seurat"))


# Visualize SEAcells annotations ----

seurat_data <- SeuratDisk::LoadH5Seurat(paste0(data_dir, "seurat/", dataset, ".with_seacells.h5Seurat"))

seacell_proportions <- 
  seurat_data@meta.data %>% 
  filter("immune" == cell_L1) %>% 
  group_by(condition, SEACell) %>% 
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(prop = n / sum(n) * 100) 

myeloid_seacells <- 
  seurat_data@meta.data %>% 
  filter("Mye" == cell_L2) %>% 
  group_by(SEACell) %>% 
  summarise(n = n(), .groups = "drop_last") %>% 
  filter(n > 200)


seacell_proportions %>%
  filter(SEACell %in% myeloid_seacells$SEACell) %>%
  ggplot(aes(SEACell, prop, fill = condition)) + 
  scale_fill_prism(palette = "colorblind_safe") + 
  geom_col(width=.64, position=position_dodge()) + 
  ggtitle(dataset) + 
  ggprism::theme_prism() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# SEACell-489 and SEACell-64

dge.seacell.489 <- 
  seurat_data %>% 
  FindMarkers(
    group.by = "SEACell",
    ident.1 = "SEACell-489",
    min.pct = 0.1,
    logfc.threshold = 0.25, 
    recorrect_umi = FALSE
  ) %>% 
  mutate(avg_log2FC = -avg_log2FC)
  rownames_to_column("gene")
  

dge.seacell.489 

dge.seacell.489 %>% filter(grepl("HSV", gene, ignore.case = TRUE)) %>% arrange(p_val_adj)

dge.seacell.64 <- 
  seurat_data %>% 
  FindMarkers(
    group.by = "SEACell",
    ident.1 = "SEACell-64",
    min.pct = 0.1,
    logfc.threshold = 0.25, 
    recorrect_umi = FALSE
  ) %>% 
  mutate(avg_log2FC = -avg_log2FC) %>%
  rownames_to_column("gene")

dge.seacell.64 

dge.seacell.64 %>% filter(grepl("HSV", gene, ignore.case = TRUE)) %>% arrange(p_val_adj)

# seurat_subset <- subset(seurat_data, subset = SEACell %in% c("SEACell-489", "SEACell-64"))
# 
# VlnPlot(seurat_subset, assay = "SCT", 
#         cols = c("black", "magenta"),
#         features = c("HSV-RS1", "HSV-UL50"), 
#         group.by = "SEACell", 
#         split.by = "condition", split.plot = TRUE,
#         stack = TRUE, flip = TRUE) +
#   ggtitle("Combined TG")
