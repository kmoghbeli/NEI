## Utility file to make GTF file from a Genbank record
library(tidyverse)

# This file takes the mouse transcripts2genes file and either the KOS or RE transcripts2genes files 
# and makes a combined transcripts2genes file for each

data_dir <- "../data_objects/"

mouse_t2g <- readr::read_tsv(paste0(data_dir, "mouse_index_standard_nov2023/t2g.txt"), col_names = FALSE)
kos_t2g <- readr::read_tsv(paste0(data_dir, "viral_genomes/kos.transcripts_to_genes.txt"), col_names = FALSE) %>% 
  select(1:3) 
re_t2g <- readr::read_tsv(paste0(data_dir, "viral_genomes/re.transcripts_to_genes.txt"), col_names = FALSE) %>% 
  select(1:3)

# Write the mouse plus kos
bind_rows(mouse_t2g, kos_t2g) %>% 
  readr::write_tsv(paste0(data_dir, "kallisto_mouse_index_v10.with_kos.transcripts_to_genes.txt"), col_names = FALSE)

# Write the mouse plus re
bind_rows(mouse_t2g, re_t2g) %>% 
  readr::write_tsv(paste0(data_dir, "kallisto_mouse_index_v10.with_re.transcripts_to_genes.txt"), col_names = FALSE)
