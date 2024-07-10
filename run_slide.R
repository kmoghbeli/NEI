#devtools::install_github("jishnu-lab/SLIDE")

library(tidyverse)
library(SLIDE)
library(doParallel)


cores <-  as.numeric(Sys.getenv('SLURM_CPUS_PER_TASK', unset=NA))

if(is.na(cores)) cores <- detectCores()

cat('number of cores using', cores, '. . .\n')

registerDoParallel(cores)

yaml_path <- commandArgs(trailingOnly=TRUE)
#yaml_path <- "./slide_runs/combined_tg_Mac_KOSvRE/er.yaml"

print(paste0("Running SLIDE on: ", yaml_path))

yaml_args <- yaml::read_yaml(yaml_path)
 
optimizeSLIDE(yaml_args, sink_file = FALSE)

#summary_table <- read.csv(paste0(yaml_args$out_path, "/summary_table.csv"), row.names = 1)

########################################################################################




## Clean up parallel stuff
registerDoSEQ()

env <- foreach:::.foreachGlobals
rm(list=ls(name=env), pos=env)

