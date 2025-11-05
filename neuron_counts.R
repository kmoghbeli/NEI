# Import Neuron Counts From Kallisto

library(tidyverse)
library(tximport)

data_dir <- "../data_objects/"

# kallisto_results <- "../data_objects/kallisto_neuron_outs/"
kallisto_results <- paste0(data_dir, "kallisto_neuron_outs.with_viral/")

metadata <- readr::read_csv(paste0(data_dir, "neuron_metadata_with_filenames.csv"))

transcript2gene.mouse <- readr::read_tsv(paste0(data_dir, "mouse_index_standard_nov2023/t2g.txt"), col_names = FALSE) %>% 
  rename_with(~ c('TXNAME', 'ENSID', 'GENEID'), 1:3) %>% 
  select(TXNAME, GENEID) %>% 
  drop_na(GENEID)

transcript2gene.mouse_with_kos <- readr::read_tsv(paste0(data_dir, "kallisto_mouse_index_v10.with_kos.transcripts_to_genes.txt"), col_names = FALSE) %>% 
  rename_with(~ c('TXNAME', 'ENSID', 'GENEID'), 1:3) %>% 
  select(TXNAME, GENEID) %>% 
  drop_na(GENEID)

transcript2gene.mouse_with_re <- readr::read_tsv(paste0(data_dir, "kallisto_mouse_index_v10.with_re.transcripts_to_genes.txt"), col_names = FALSE) %>% 
  rename_with(~ c('TXNAME', 'ENSID', 'GENEID'), 1:3) %>% 
  select(TXNAME, GENEID) %>% 
  drop_na(GENEID)



files <- list.files(path = kallisto_results, pattern = "abundance.tsv", recursive = TRUE)

filenames_metadata <- 
  files %>% as_tibble() %>% 
  rename(file = value) %>% 
  mutate(cell = sub("/abundance.tsv", "", file), 
         name = paste0("neuron_", cell),
         file = normalizePath(paste0(kallisto_results, file))) %>% 
  left_join(metadata %>% select(cell, condition), 
            by = join_by("cell" == "cell")) %>% 
  select(name, cell, condition, file)

# Ref: https://support.bioconductor.org/p/132550/
txi.control_and_scratch <- 
  tximport(filenames_metadata %>% filter(condition %in% c("control", "scratch")) %>% arrange(cell) %>% pull(file), 
           type = "kallisto", 
           tx2gene = transcript2gene.mouse, 
           ignoreAfterBar = TRUE,
           countsFromAbundance = "lengthScaledTPM")

txi.kos <- 
  tximport(filenames_metadata %>% filter(condition %in% c("kos")) %>% arrange(cell) %>% pull(file), 
           type = "kallisto", 
           tx2gene = transcript2gene.mouse_with_kos, 
           ignoreAfterBar = TRUE,
           countsFromAbundance = "lengthScaledTPM")

txi.re <- 
  tximport(filenames_metadata %>% filter(condition %in% c("re")) %>% arrange(cell) %>% pull(file), 
           type = "kallisto", 
           tx2gene = transcript2gene.mouse_with_re, 
           ignoreAfterBar = TRUE,
           countsFromAbundance = "lengthScaledTPM")

neuron_counts.control_and_scratch <- 
  txi.control_and_scratch$counts %>% 
  t() %>% 
  as_tibble() %>% 
  mutate(cell = filenames_metadata %>% filter(condition %in% c("control", "scratch")) %>% arrange(cell) %>% pull(name)) %>% 
  column_to_rownames("cell")

neuron_counts.kos <- 
  txi.kos$counts %>% 
  t() %>% 
  as_tibble() %>% 
  mutate(cell = filenames_metadata %>% filter(condition %in% c("kos")) %>% arrange(cell) %>% pull(name)) %>% 
  column_to_rownames("cell")

neuron_counts.re <- 
  txi.re$counts %>% 
  t() %>% 
  as_tibble() %>% 
  mutate(cell = filenames_metadata %>% filter(condition %in% c("re")) %>% arrange(cell) %>% pull(name)) %>% 
  column_to_rownames("cell")

all_counts <- 
  bind_rows(neuron_counts.control_and_scratch,
            neuron_counts.kos, 
            neuron_counts.re)

# Set any NA values to 0
all_counts[is.na(all_counts)] <- 0

# Round the non-integer counts
all_counts <- round(all_counts)

# Transpose the data frame
all_counts <- all_counts %>% 
  t() %>% 
  as_tibble(rownames = "gene")

all_counts %>% readr::write_csv("../../../_DATASETS/NeuroImmune/neuron_counts_all_with_viral_kallisto_tximport_lengthscaledtpm.csv")


 
################################################################################################
################################################################################################
################################################################################################

# neuron_counts <- readr::read_csv("../../../_DATASETS/NeuroImmune/control_corneal_afferent_counts.csv") %>% 
#   rename(ENSG = `...1`) %>% 
#   mutate(ENSG = sub("\\.\\d+$", "", ENSG))   # remove Ensembl Version number from the Ensembl Gene ID
# 
# ## Ensembl-centric Database to map Ensembl IDs to Gene symbols
# #library(EnsDb.Mmusculus.v79)
# 
# ensg_to_gene <- AnnotationDbi::select(EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79,
#                                       keys = neuron_counts$ENSG,
#                                       keytype = "GENEID",
#                                       columns = c("SYMBOL"))
# 
# neuron_counts_filtered <- neuron_counts %>% 
#   left_join(ensg_to_gene, by = join_by(ENSG == GENEID)) %>% 
#   dplyr::filter(!is.na(SYMBOL) & "" != SYMBOL) %>% 
#   dplyr::rename(gene = SYMBOL) %>% 
#   dplyr::select(-c(ENSG)) %>% 
#   group_by(gene) %>% 
#   summarise(across(everything(), sum)) %>% 
#   select(gene, where(~ is.numeric(.) && sum(.) > 0))   # Don't keep any cells/neurons with no counts for any genes
# 
# # Round the non-integer counts
# neuron_counts_filtered <- round(neuron_counts_filtered %>% column_to_rownames("gene")) %>% rownames_to_column("gene")
# 
# neuron_counts_filtered %>% readr::write_csv("../../../_DATASETS/NeuroImmune/control_corneal_afferent_counts_filtered.csv")
# 
# neuron_counts_filtered %>% 
#   pivot_longer(-gene) %>%   
#   pivot_wider(names_from = gene, values_from = value) %>% 
#   readr::write_csv("../../../_DATASETS/NeuroImmune/control_corneal_afferent_counts_filtered_transposed.csv")


