#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#### SCRIPT INTRODUCTION ####
#
#' @name _targets.R  
#' @description R script to launch the target pipeline
#' @author Anne Baranger
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Options and packages ----------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Load targets
library(targets)
# library(tarchetypes)
# Load functions
# lapply(grep("R$", list.files("R"), value = TRUE), function(x) source(file.path("R", x)))
source("R/functions_data.R")
source("R/functions_analyses.R")


# install if needed and load packages
packages.in <- c("stringr","data.table","tidyr","future",
                 "terra",
                 "dplyr","sf","rstan")
#for(i in 1:length(packages.in)) if(!(packages.in[i] %in% rownames(installed.packages()))) install.packages(packages.in[i])
# Targets options
options(tidyverse.quiet = TRUE, clustermq.scheduler = "multiprocess")
tar_option_set(packages = packages.in,
               memory = "transient")
future::plan(future::multisession, workers = 6)
set.seed(2)
occurence = readRDS("occurence.binom")
df.species= readRDS("df.species")
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Targets workflow --------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

list(
  tar_target(
    species_fit_list,
    df.species |> filter(!is.na(px),!is.na(lt50)) |> 
      filter(species %in% unique(occurence$species)) |> 
      pull(species)
  ),
  tar_target(
    mod_random,
    fit_random_sp(occurence,
                  df.species,
                  var.hsm="psi_eraday_real",
                  var.fsm="tmin_era",
                  sp.excl=species_fit_list,
                  folder.out="mod.rdata"),
    pattern=map(species_fit_list),
    iteration="vector"
  ),
  NULL
)


