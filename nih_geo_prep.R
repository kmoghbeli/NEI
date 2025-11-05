library(tidyverse)

base_dir <- "/ix/cigcore/shared/cigcore_wam30_rfe4/"

# Make a filename dump for NIH GEO upload of all the Neuron FASTQ files (R1 and R2) ----
neuron_metadata <- readr::read_csv("../data_objects/neuron_metadata_with_filenames.csv")
neuron_metadata %>% filter(condition == "control") %>% select(filename) %>% 
  mutate(filename2 = sub("_R1_", "_R2_", filename)) %>% 
  pivot_longer(cols = c(filename, filename2), values_to = "filename") %>% 
  select(-name) %>% 
  arrange(filename) %>%
  pull(filename) %>% 
  paste(collapse = "\n") %>% 
  readr::write_file("../../Manuscripts/Cell Reports - Corneal Afferents/GEO Submission Data/neuron_filenames_for_GEO.txt")

# Extract just the control neurons from the neuron counts and write them to disk
counts_filename <- "../../../_DATASETS/NeuroImmune/neuron_counts_all_kallisto_tximport_lengthscaledtpm.csv"

control_neuron_names <- 
  neuron_metadata %>% 
  filter(condition == "control") %>% 
  mutate(cell = paste0("neuron_", cell)) %>% 
  pull(cell)

readr::read_csv(counts_filename, show_col_types = FALSE) %>% 
  column_to_rownames(var = "gene") %>% 
  select(all_of(control_neuron_names)) %>% 
  rownames_to_column(var = "gene") #%>%
  readr::write_csv("../../Manuscripts/Cell Reports - Corneal Afferents/GEO Submission Data/control_neuron_counts_kallisto_tximport_lengthscaledtpm.csv")
  
