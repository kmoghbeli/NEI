library(tidyverse)

base_dir <- "/ix/cigcore/shared/cigcore_wam30_rfe4/"

get_filename <- function(base_dir, dir, prefix, num) {
  
  filename <- c()
  
  filename <- list.files(path = paste0(base_dir, dir), 
                         pattern = paste0("^", prefix, num, "_.+R1_001.fastq.gz"), 
                         full.names = FALSE)
  
  if (0 == length(filename)) {
    filename <- list.files(path = paste0(base_dir, dir), 
                           pattern = paste0("^", prefix, "_", num, "_.+R1_001.fastq.gz"), 
                           full.names = FALSE)
  }
  
  if (0 == length(filename)) {
    warning("No matching file found.")  
    filename <- ""
  } else if (length(filename) > 1) {
    warning("More than 1 file matched when only expecting 1 - keeping first. ", filename)
    filename <- filename[1]
  } 
  
  return(filename)
}

### FILENAME PREP - ONLY NEED TO DO THIS ONCE
## Read in Excel file with neuron naming info
neuron_filenames <- readr::read_csv("data/neuron_metadata.csv") %>% 
  mutate(cell_tube_prefix = ifelse(is.na(cell_tube_prefix), "", cell_tube_prefix), 
         cell_tube_begin = as.numeric(str_extract(cell_tube_range, "^\\d+")), 
         cell_tube_end = as.numeric(str_extract(cell_tube_range, "\\d+$")), 
         n = cell_tube_end - cell_tube_begin + 1) %>% 
  uncount(n, .id = "rel_id") %>% 
  rowwise() %>% 
  mutate(cell_tube_num = cell_tube_begin + rel_id - 1, 
         cell = paste0(cell_tube_prefix, cell_tube_num), 
         filename = paste0(base_dir, folder, "/", get_filename(base_dir, folder, cell_tube_prefix, cell_tube_num))) %>% 
  rename(batch_folder = folder) %>% 
  select(condition, cell_tube, cell, date, species, animal_num, sex, batch_folder, filename) 

## Write a new CSV file with filenames for each individual neuron
neuron_filenames %>% 
  readr::write_csv("neuron_metadata_with_filenames.csv")

neuron_filenames %>% pull(filename) %>% paste(collapse = " ") %>% readr::write_file("neuron_filenames_only.txt")

## Create Kallisto SLURM file
kallisto_template <- readr::read_lines("neuron_kallisto.slurm.template")
kallisto_commands <- readr::read_csv("neuron_metadata_with_filenames.csv") %>% 
  mutate(kallisto_command = paste0("kallisto quant -i /ix/djishnu/kaveh/kallisto.mouse_index.v10/transcriptome.idx -o kallisto_neuron_outs/", 
                                   cell, " ", filename, " ", sub("_R1_", "_R2_", filename))) %>% 
  pull(kallisto_command)

readr::write_lines(c(kallisto_template, "\n", kallisto_commands), "neuron_kallisto.slurm")


