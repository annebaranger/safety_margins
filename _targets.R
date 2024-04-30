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
source("R/functions_data.R")
source("R/functions_analyses.R")


# install if needed and load packages
packages.in <- c("stringr","ggplot2","data.table","tidyr",
                 "viridis","rgdal","raster","rosm","terra","rnaturalearth",
                 "dplyr","sf","lubridate","lme4")
for(i in 1:length(packages.in)) if(!(packages.in[i] %in% rownames(installed.packages()))) install.packages(packages.in[i])
# Targets options
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
  
  # Soil water content from ERA5, monthly timestep, 4 horizons, perc05
  tar_target(
    swc_era_month,
    get_swc_grib(dir.data="data/ERA5-land/monthly/",
                 dir.file="era5_svwcperc05_",
                 vars=c("h1","h2","h3","h4"),
                 extension="_1950-2021_inv.grib",
                 europe)
  ),
  
    # Weight horizons until 289cm
    tar_target(
      swc_month_289,
      weight_swc(swc_era_month,
                 depth=289)
    ),
  
    # Weight horizons until 100cm
    tar_target(
      swc_month_100,
      weight_swc(swc_era_month,
                 depth=100)
    ),
  
  # Soil water content from ERA5, daily timestep, 4 horizons, mean of 5 min
  tar_target(
    swc_era_day,
    get_swc_csv(dir.data="data/ERA5-land/daily/",
                dir.file="swcd-1950-2021-",
                vars=c("layer1","layer2","layer3","layer4"),
                extension="_min.csv",
                rast.model="data/ERA5-land/daily/era5_svwcperc05d_h1_1984-2021_inv.grib",
                europe)
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
  
  # soil water content from CERRA whole year, daily timestep, 11 horizons, perc05
  tar_target(
    swc_cerra_day,
    get_swc_grib(dir.data="data/cerra-land/liquid_vol_content/min/",
                 dir.file="cerra_lvslperc05_",
                 vars=paste0("h",1:11),
                 extension="_1984-2021_inv.grib",
                 europe)
  ),
  
  # soil water content from CERRA summer only, daily timestep, 11 horizons, perc05
  tar_target(
    swc_cerra_day_summer,
    get_swc_grib(dir.data="data/cerra-land/liquid_vol_content/min_sum/",
                 dir.file="cerra_lvslperc05_",
                 vars=paste0("h",1:11),
                 extension="_1984-2021_inv.grib",
                 europe)
  ),
  
  # min temperature from cerra, skin temp

  tar_target(
    tmin_cerra,
    get_frostindex_cerra(europe,
                         dir.file="data/cerra-land/skin_temp/cerra_perc05_1984-2021_inv.grib")
  ),
  
  # min temperature from cerra, 2m
  tar_target(
    tmin2m_cerra,
    get_frostindex_cerra(europe,
                         dir.file="data/cerra-land/skin_temp/cerra_t2mperc05_1984-2021_inv.grib")
  ),
  
  # min temprature from ERA5
  tar_target(
    tmin2m_era,
    {dir.file="data/ERA5-land/t2m/era5_t2min05_1984-2021.grib"
    tmin=rast(dir.file)-273.15
    return(as.data.frame(tmin,xy=TRUE))
    }
  ),

  # min temprature from chelsa
  tar_target(
    tmin_chelsa,
    get_frostindex_chelsa(europe,
                          dir.file="data/CHELSA/CHELSA_EUR11_tasmin_month_min_19802005.nc")
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
    psi_eramonth_real,
    compute_psi_sureau(swc_era_month,
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
                       "output/psi_era_month_real_fixed.csv"
    )
  ),
  tar_target(
    psi_eraday_real,
    compute_psi_sureau(swc_era_day,
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
                       "output/psi_era_day_real_min5_fixed.csv"
    )
  ),
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
  tar_target(
    psi_eraday_100,
    compute_psi_sureau(swc_era_day,
                       europe,
                       dir.hydro="data/EU_SoilHydroGrids_1km/",
                       depth_max=100,
                       dir.depth="data/STU_EU_Layers/STU_EU_DEPTH_ROOTS.rst",
                       dir.ecoregions="data/WWF/official",
                       LAImax=5,
                       fRootToLeaf=1,
                       rootRadius=0.0004,
                       beta=0.97,
                       obs=c(0,0.07,0.28,1,2.89),
                       ref=c(0,0.05,0.15,0.3,0.6,1,2),
                       max_depth=3,
                       "output/psi_era_day_100_fixed.csv"
    )
  ),
  tar_target(
    psi_cerraday_real,
    compute_psi_sureau(swc_cerra_day,
                       europe,
                       dir.hydro="data/EU_SoilHydroGrids_1km/",
                       depth_max=NULL,
                       dir.depth="data/STU_EU_Layers/STU_EU_DEPTH_ROOTS.rst",
                       dir.ecoregions="data/WWF/official",
                       LAImax=5,
                       fRootToLeaf=1,
                       rootRadius=0.0004,
                       beta=0.97,
                       obs=c(0,0.01,0.04,0.1,0.2,0.4,0.6,0.8,1,1.5,2,3),
                       ref=c(0,0.05,0.15,0.3,0.6,1,2),
                       max_depth=3,
                       "output/psi_cerra_day_real_fixed.csv"
    )
  ),
  tar_target(
    psi_cerraday_100,
    compute_psi_sureau(swc_cerra_day,
                       europe,
                       dir.hydro="data/EU_SoilHydroGrids_1km/",
                       depth_max=100,
                       dir.depth="data/STU_EU_Layers/STU_EU_DEPTH_ROOTS.rst",
                       dir.ecoregions="data/WWF/official",
                       LAImax=5,
                       fRootToLeaf=1,
                       rootRadius=0.0004,
                       beta=0.97,
                       obs=c(0,0.01,0.04,0.1,0.2,0.4,0.6,0.8,1,1.5,2,3),
                       ref=c(0,0.05,0.15,0.3,0.6,1,2),
                       max_depth=3,
                       "output/psi_cerra_day_100_fixed.csv"
    )
  ),
  # tar_target(
  #   psi_cerraday_real_beta,
  #   compute_psi_sureau(swc_cerra_day,
  #                      europe,
  #                      dir.hydro="data/EU_SoilHydroGrids_1km/",
  #                      depth_max=NULL,
  #                      dir.depth="data/STU_EU_Layers/STU_EU_DEPTH_ROOTS.rst",
  #                      dir.ecoregions="data/WWF/official",
  #                      LAImax=5,
  #                      fRootToLeaf=1,
  #                      rootRadius=0.0004,
  #                      beta=NULL,
  #                      obs=c(0,0.01,0.04,0.1,0.2,0.4,0.6,0.8,1,1.5,2,3),
  #                      ref=c(0,0.05,0.15,0.3,0.6,1,2),
  #                      max_depth=3,
  #                      "output/psi_cerra_day_real_beta.csv"
  #   )
  # ),
  tar_target(
    psi_cerradaysum_real,
    compute_psi_sureau(swc_cerra_day_summer,
                       europe,
                       dir.hydro="data/EU_SoilHydroGrids_1km/",
                       depth_max=NULL,
                       dir.depth="data/STU_EU_Layers/STU_EU_DEPTH_ROOTS.rst",
                       dir.ecoregions="data/WWF/official",
                       LAImax=5,
                       fRootToLeaf=1,
                       rootRadius=0.0004,
                       beta=0.97,
                       obs=c(0,0.01,0.04,0.1,0.2,0.4,0.6,0.8,1,1.5,2,3),
                       ref=c(0,0.05,0.15,0.3,0.6,1,2),
                       max_depth=3,
                       "output/psi_cerra_day_sum_real_fixed.csv"
    )
  ),
  # tar_target(
  #   psi_cerradaysum_100,
  #   compute_psi_sureau(swc_cerra_day_summer,
  #                      europe,
  #                      dir.hydro="data/EU_SoilHydroGrids_1km/",
  #                      depth_max=100,
  #                      dir.depth="data/STU_EU_Layers/STU_EU_DEPTH_ROOTS.rst",
  #                      dir.ecoregions="data/WWF/official",
  #                      LAImax=5,
  #                      fRootToLeaf=1,
  #                      rootRadius=0.0004,
  #                      beta=0.97,
  #                      obs=c(0,0.01,0.04,0.1,0.2,0.4,0.6,0.8,1,1.5,2,3),
  #                      ref=c(0,0.05,0.15,0.3,0.6,1,2),
  #                      max_depth=3,
  #                      "output/psi_cerra_day_100_fixed.csv"
  #   )
  # ),
  # tar_target(
  #   psi_cerradaysum_real_beta,
  #   compute_psi_sureau(swc_cerra_day_summer,
  #                      europe,
  #                      dir.hydro="data/EU_SoilHydroGrids_1km/",
  #                      depth_max=NULL,
  #                      dir.depth="data/STU_EU_Layers/STU_EU_DEPTH_ROOTS.rst",
  #                      dir.ecoregions="data/WWF/official",
  #                      LAImax=5,
  #                      fRootToLeaf=1,
  #                      rootRadius=0.0004,
  #                      beta=NULL,
  #                      obs=c(0,0.01,0.04,0.1,0.2,0.4,0.6,0.8,1,1.5,2,3),
  #                      ref=c(0,0.05,0.15,0.3,0.6,1,2),
  #                      max_depth=3,
  #                      "output/psi_cerra_day_real_beta.csv"
  #   )
  # ),
  
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
    df.LT50,
    get_LT50(file.output="output/df_LT50_filtered.csv")
  ),
  tar_target(
    df.P50,
    get_P50(species.list,
            file.output="output/df_P50_filtered.csv")
  ),
  tar_target(
    df.traits,
    get_traits(df.P50=df.P50,
               df.LT50=df.LT50$df.LT50sp.cor,
               file.output="output/df_trait_filtered.csv")
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
  NULL
)


