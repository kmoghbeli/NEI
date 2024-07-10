library(tidyverse)
library(EssReg)
library(doParallel)


cores <-  as.numeric(Sys.getenv('SLURM_CPUS_PER_TASK', unset=NA))

if(is.na(cores)) cores <- detectCores()

cat('number of cores using', cores, '. . .\n')

registerDoParallel(cores)

yaml_path <- commandArgs(trailingOnly=TRUE)
#yaml_path <- "./slide_runs/combined_tg_Mac_KOSvRE_D0.01_L1.0_S0.5/er.yaml"

print(paste0("Running pipeline 3 on: ", yaml_path))
 
pipelineER3(yaml_path)

yaml_args <- yaml::yaml.load_file(yaml_path)

# Set a default "spec"
yaml_args$spec <- ifelse(is.null(yaml_args$spec), 0.1, yaml_args$spec)

# check for ER results in output folder
er_results_path = list.files(yaml_args$out_path, pattern = "final_delta")

if (length(er_results_path) > 0) {
  print("ER completed successfully.")
} else {
  print("ER did not complete successfully.")
  return(-1)
}

## Clean up parallel stuff
registerDoSEQ()

env <- foreach:::.foreachGlobals
rm(list=ls(name=env), pos=env)

