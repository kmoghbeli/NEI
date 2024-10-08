# Fix Neuron Counts

library(tidyverse)
library(tximport)

kallisto_results <- "kallisto_neuron_outs"


# Ensembl Transcriptomes v96
transcript2gene <- readr::read_tsv("transcripts_to_genes.txt", col_names = FALSE) %>% 
  rename_with(~ c('TXNAME', 'ENSID', 'GENEID'), 1:3) %>% 
  select(TXNAME, GENEID)

files <- list.files(path = "kallisto_neuron_outs", pattern = "*.h5", recursive = TRUE)
cellnames <- sub("/abundance.h5", "", files)
files <- normalizePath(paste0("kallisto_neuron_outs/", files))
names(files) <- paste0("neuron_", cellnames)

# Ref: https://support.bioconductor.org/p/132550/
txi <- tximport(files, type = "kallisto", tx2gene = transcript2gene, countsFromAbundance = "lengthScaledTPM")

neuron_counts <- txi$counts %>% as_tibble(rownames = "gene") %>% 
  select(str_sort(colnames(neuron_counts), numeric = TRUE)) %>% 
  column_to_rownames("gene") %>% 
  round() %>%   # round the non-integer counts 
  rownames_to_column("gene")

neuron_counts %>% readr::write_csv("neuron_counts_all_kallisto_tximport_lengthscaledtpm.csv")

 
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


