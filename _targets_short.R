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
library(tarchetypes)
# Load functions
# lapply(grep("R$", list.files("R"), value = TRUE), function(x) source(file.path("R", x)))
source("functions_data.R")
source("functions_analyses.R")


# install if needed and load packages
packages.in <- c("stringr","ggplot2","data.table","tidyr","future",
                 "viridis","raster","rosm","terra","rnaturalearth", #"rgdal",
                 "dplyr","sf","lubridate","lme4","rstan","pROC")
options(tidyverse.quiet = TRUE, clustermq.scheduler = "multiprocess")
tar_option_set(packages = packages.in,
               memory = "transient")
future::plan(future::multisession, workers = 6)
set.seed(2)

LAImax=c(2,3,4,5,5,5,5,6,7,8)
beta=c(0.966,0.966,0.966,0.914,0.942,0.966,0.976,0.966,0.966,0.966)
values <- data.frame(LAImax=LAImax,
                     beta=beta,
                     output=paste0("output/psi_era_day_real_",LAImax,"_",beta))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Targets workflow --------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

list(
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Load data -----
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  # Europe extent
  tar_target(
    europe,
    get_europe()
  ),
 # soil water content from ERA5, daily timestep, 4 horizons, perc05
  tar_target(
    swc_era_day_05,
    get_swc_grib(dir.data="data/ERA5-land/daily/",
                 dir.file="era5_svwcperc05d_",
                 vars=c("h1","h2","h3","h4"),
                 extension="_1984-2021_inv.grib",
                 europe)
  ),
  
  # min temprature from ERA5
  tar_target(
    tmin2m_era,
    {dir.file="data/ERA5-land/t2m/era5_t2min05_1984-2021.grib"
    tmin=rast(dir.file)-273.15
    return(as.data.frame(tmin,xy=TRUE))
    }
  ),
  # extract chelsa data and scale it to era resolution
  tar_target(
    clim_chelsa,
    get_waisgdd(dir.chelsa="data/CHELSA/",
                file=c("bio1","bio12","pet_penman_mean","gdd5"),
                rast.mod="data/ERA5-land/t2m/era5_t2min05_1984-2021.grib")
  ),
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Compute psi_min -----
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  tar_target(
    psi_eradaycdo_real,
    compute_psi_sureau(swc_era_day_05,
                       europe,
                       dir.hydro="data/EU_SoilHydroGrids_1km/",
                       depth_max=NULL,
                       dir.depth="data/STU_EU_Layers/STU_EU_DEPTH_ROOTS.rst",
                       dir.ecoregions="data/WWF/official",
                       LAImax=5,
                       fRootToLeaf=1,
                       rootRadius=0.0004,
                       beta=0.97,
                       obs=c(0,0.07,0.28,1,2.89),
                       ref=c(0,0.05,0.15,0.3,0.6,1,2),
                       max_depth=3,
                       "output/psi_era_day_real_fixed.csv"
    )
  ),
  tar_target(
    psi_eradaycdo_real_beta,
    compute_psi_sureau(swc_era_day_05,
                       europe,
                       dir.hydro="data/EU_SoilHydroGrids_1km/",
                       depth_max=NULL,
                       dir.depth="data/STU_EU_Layers/STU_EU_DEPTH_ROOTS.rst",
                       dir.ecoregions="data/WWF/official",
                       LAImax=5,
                       fRootToLeaf=1,
                       rootRadius=0.0004,
                       beta=NULL,
                       obs=c(0,0.07,0.28,1,2.89),
                       ref=c(0,0.05,0.15,0.3,0.6,1,2),
                       max_depth=3,
                       "output/psi_era_day_real_beta.csv"
    )
  ),
  
  # sensitivity analysis
  tar_map(
    values=values,
    tar_target(sensitivity,
               compute_psi_sureau(swc_era_day_05,
                                  europe,
                                  dir.hydro="data/EU_SoilHydroGrids_1km/",
                                  depth_max=NULL,
                                  dir.depth="data/STU_EU_Layers/STU_EU_DEPTH_ROOTS.rst",
                                  dir.ecoregions="data/WWF/official",
                                  LAImax,
                                  fRootToLeaf=1,
                                  rootRadius=0.0004,
                                  beta,
                                  obs=c(0,0.01,0.04,0.1,0.2,0.4,0.6,0.8,1,1.5,2,3),
                                  ref=c(0,0.05,0.15,0.3,0.6,1,2),
                                  max_depth=3,
                                  output
               )
    )
  ),
  tar_target(
    clim_chelsa_node,
    clim_chelsa |> 
      mutate(cell=row_number())
  ),
  tar_target(
    psi_eradaycdo_real_dws,
    terra::project(rast(psi_eradaycdo_real[,c("x","y","psi")],
                        crs="epsg:4326"),
                   rast(clim_chelsa,crs="epsg:4326"),
                   method="average") |> 
      as.data.frame(xy=TRUE)
  ),
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Load traits -----
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  tar_target(
    species.list,
    list.files("data/chorological_maps_dataset/")
  ),
  tar_target(
    traits_max,
    get_traits_max(species.list,
                   dir.file="data/Species traits/base_traits_P50_LTx.xlsx",
                   file.output="output/traits_filtered_max.csv")
  ),
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Build occurence dataset -----
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  tar_target(
    db.mauri,
    get_mauri(species.list,
              dir.occ="data/EUForestsMauri/EUForestspecies.csv")
  ),
  
  tar_target(
    margin.limit,
    get.margins(occurence,
                psi="psi_eraday_real",
                tmin="tmin_era")
  ),
  
  
  tar_target(
    occurence,
    get.occ.clim(db.mauri,
                 clim.list=list(psi_eraday_real=psi_eradaycdo_real[,c("x","y","psi")],
                                tmin_era=tmin2m_era,
                                mat=clim_chelsa[,c("x","y","mat")],
                                map=clim_chelsa[,c("x","y","map")],
                                pet=clim_chelsa[,c("x","y","pet")],
                                sgdd=clim_chelsa[,c("x","y","sgdd")],
                                wai=clim_chelsa[,c("x","y","wai")]),
                 dws=FALSE,
                 file.path="output/df.occurence.clim.csv")
  ),
  tar_target(
    occurence.binom,
    get.occ.clim(db.mauri,
                 clim.list=list(psi_eraday_real=psi_eradaycdo_real_dws[,c("x","y","psi")],
                                tmin_era=tmin2m_era,
                                mat=clim_chelsa_node[,c("x","y","mat")],
                                map=clim_chelsa_node[,c("x","y","map")],
                                pet=clim_chelsa_node[,c("x","y","pet")],
                                sgdd=clim_chelsa_node[,c("x","y","sgdd")],
                                wai=clim_chelsa_node[,c("x","y","wai")],
                                cell=clim_chelsa_node[,c("x","y","cell")]),
                 dws=TRUE,
                 file.path="output/df.occurence.clim.binom.csv")
  ),
  tar_target(
    occurence_beta,
    get.occ.clim(db.mauri,
                 clim.list=list(psi_eraday_real=psi_eradaycdo_real[,c("x","y","psi")],
                                psi_eraday_real_beta=psi_eradaycdo_real_beta[,c("x","y","psi")],
                                tmin_era=tmin2m_era,
                                mat=clim_chelsa[,c("x","y","mat")],
                                map=clim_chelsa[,c("x","y","map")],
                                pet=clim_chelsa[,c("x","y","pet")],
                                sgdd=clim_chelsa[,c("x","y","sgdd")],
                                wai=clim_chelsa[,c("x","y","wai")]),
                 dws=FALSE,
                 file.path="output/df.occurence.clim.beta.csv")
  ),
  tar_target(
    df.preval,
    get.prevalence(species.list,
                   occurence)
  ),
  tar_target(
    df.niche,
    get.niche(species.list,
              occurence,
              psi="psi_eraday_real",
              tmin="tmin_era")
  ),
  tar_target(
    df.shadetol,
    get.shadetol(species.list,
                 occurence)
  ),
  tar_target(
    df.species,
    get.species(species.list,
                df.preval,
                df.shadetol,
                df.niche,
                traits_max,
                file.output="output/df.species2.csv")
  ),
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Fit species model  -----
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  tar_target(
    species_fit_list,
    df.species |> filter(!is.na(px),!is.na(lt50)) |> 
      filter(species %in% unique(occurence$species)) |> 
      pull(species)
  ),
  tar_target(
    fit.safmarg.sp,
    fit_species(occurence=occurence,
                var.hsm="psi_eraday_real",
                var.fsm="tmin_era",
                df.species=df.species,
                sp=species_fit_list,
                output="fit_safmarg_era/"),
    pattern=map(species_fit_list)
  ),
  
  tar_target(
    fit.safmarg.sp.beta,
    fit_species(occurence=occurence_beta,
                var.hsm="psi_eraday_real",
                var.fsm="tmin_era",
                df.species=df.species,
                sp=species_fit_list,
                output="fit_safmarg_era_beta/"),
    pattern=map(species_fit_list)
  ),
  tar_target(
    safmarg_sp_select,
    fit.safmarg.sp %>% 
      filter(rhat<1.2) %>% 
      filter(divergence <0.1) %>% 
      group_by(species) %>% 
      slice(which.min(bic)) %>% 
      ungroup()
  ),
  tar_target(
    safmarg_spbeta_select,
    fit.safmarg.sp.beta %>% 
      filter(rhat<1.2) %>% 
      filter(divergence <0.1) %>% 
      group_by(species) %>% 
      slice(which.min(bic)) %>% 
      ungroup()
  ),
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Fit generic model  -----
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  tar_target(
    mod_random,
    fit_random_sp(occurence.binom,
                  df.species,
                  var.hsm="psi_eraday_real",
                  var.fsm="tmin_era",
                  sp.excl=species_fit_list,
                  folder.out="mod.rdata/final"),
    pattern=map(species_fit_list),
    iteration="vector"
  ),
  tar_target(
    mod_random_predict,
    predict_fit(occurence.binom,
                df.species,
                var.hsm="psi_eraday_real",
                var.fsm="tmin_era",
                folder.out="mod.rdata/final")
  ),
  NULL
)


