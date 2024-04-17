#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#### SCRIPT INTRODUCTION ####
#
#' @name functions_data.R  
#' @description R script containing all functions relative to data
#               importation, subsetting and computation of 
#               required variables
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



#' Europe extent
#' 
#' @description This is a function to extract countries borders from an R 
#' package 
#' @param countries default selection encompass geographical definition of 
#' europe, except Russia
#' @return SpatialPolygonsDataFrame containing all countries of selection
#'
get_europe <- function(){
  # Get world countries data
  world <- ne_countries(scale = "medium", returnclass = "sf")
  europe <- world[world$continent == "Europe", ]
  europe <-  crop(vect(europe),ext(-10,40,32,72))
  europe <- st_as_sf(europe)
  return(europe)
}



#' Load min computed externally with python, from ERA5land
#' 
#' @description Load and reshape SWC data computed with python program, per horizon
#' @note Python code need to be run prior executing this function
#' @param dir.data
#' @param dir.file name of file
#' @param vars list of vars corresponding to horizons
#' @param extension name of the file extension
#' @param rast.model directory of a raster model
#' @param europe europe shape for croping
#' @return dataframe with min swc per horizon

get_swc_csv <- function(dir.data="data/ERA5-land/daily/",
                        dir.file="swcd-1950-2021-",
                        vars=c("layer1","layer2","layer3","layer4"),
                        extension="_min.csv",
                        rast.model="data/ERA5-land/daily/era5_svwcperc05d_h1_1984-2021_inv.grib",
                        europe){
  rast.model=rast(rast.model)
  rast.swc=rast(nrows=dim(rast.model)[1],
                ncol=dim(rast.model)[2],
                xmin=ext(rast.model)$xmin,
                xmax=ext(rast.model)$xmax,
                ymin=ext(rast.model)$ymin,
                ymax=ext(rast.model)$ymax,
                nlyrs=4)
    for (i in 1:length(vars)){
      rast.min=as.matrix(read.table(paste0(dir.data,dir.file,vars[i],extension),
                                    header=FALSE,
                                    sep=",",
                                    dec = "."))
      rast.min[rast.min==-32767] <- NA
      rast.swc[[i]]=rast(nrows=dim(rast.min)[1],
                         ncol=dim(rast.min)[2],
                         xmin=ext(rast.model)$xmin,
                         xmax=ext(rast.model)$xmax,
                         ymin=ext(rast.model)$ymin,
                         ymax=ext(rast.model)$ymax,
                         crs=crs(rast.model),
                         vals=c(t(rast.min)))
      
    }
  names(rast.swc)=paste0("h",1:length(vars))
  rast.swc=crop(project(rast.swc,
                       "EPSG:4326"),
               vect(europe))
  return(as.data.frame(rast.swc,xy=TRUE))
}



#' Load min computed externally with python, from ERA5land
#' 
#' @description Load and reshape SWC data computed with python program, per horizon
#' @note Python code need to be run prior executing this function
#' @param dir.data
#' @param dir.file name of file
#' @param vars list of vars corresponding to horizons
#' @param extension name of the file extension
#' @param rast.model directory of a raster model
#' @param europe europe shape for croping
#' @return dataframe with min swc per horizon

get_swc_grib <- function(dir.data,
                        dir.file,
                        vars,
                        extension,
                        europe){
  rast.swc=rast(paste0(dir.data,dir.file,vars,extension))
  names(rast.swc)=paste0("h",1:length(vars))
  rast.swc=crop(project(rast.swc,
                        "EPSG:4326"),
                vect(europe))
  return(as.data.frame(rast.swc,xy=TRUE))
}


#' Load frost index from cerra
#' 
#' @description Load files of minimum temperature computed with cdo
#' @param europe 
#' @param dir.file
#' @return dataframe of minimum temperature
get_frostindex_cerra<-function(europe,
                               dir.file){
  tmin<-rast(dir.file)
  tmin<-crop(project(tmin,
                     "epsg:4326"),
             vect(europe))
  names(tmin)="tmin"
  return(as.data.frame(tmin,xy=TRUE))
  
}


#' Load frost index from chelsa
#' 
#' @description Load files of minimum temperature computed with cdo
#' @param europe 
#' @param dir.file
#' @return dataframe of minimum temperature
get_frostindex_chelsa<-function(europe,
                                dir.file="data/CHELSA/CHELSA_EUR11_tasmin_month_min_19802005.nc"){
  tmin=min(rast("data/CHELSA/CHELSA_EUR11_tasmin_month_min_19802005.nc"),na.rm=FALSE)
  tmin=classify(tmin, cbind(6553500, NA)) #set as NA default value
  names(tmin)="tmin"
  return(as.data.frame(tmin,xy=TRUE))
  
}


#' Extract from Chelsea and compute WAI/SGDD
#' 
#' @description extract required variables from chelsea and computed sgdd and wai
#' on points entered as argument
#' @note Chelsea dwonloaded from website
#' @param dir.chelsa
#' @param df.loc data.frame of x and y location
#' @return df.loc + extracted variables + computed variables

get_waisgdd <- function(dir.chelsa="data/CHELSA/",
                        file=c("bio1","bio12","pet_penman_mean","gdd5"),
                        rast.mod){
  rast.mod=project(rast(rast.mod),
                   "EPSG:4326")
  clim.files=paste0(dir.chelsa,"CHELSA_",file,"_1981-2010_V.2.1.tif")
  rast.clim=rast(clim.files)
  names(rast.clim)=file
  rast.clim=resample(crop(rast.clim,
                          rast.mod),
                     rast.mod,
                     method="near")
  names(rast.clim)=c("mat","map","pet","sgdd")
  rast.clim=as.data.frame(rast.clim,xy=TRUE) |> 
    mutate(wai=(map-12*pet)/pet)
  return(rast.clim)
}


#' Compute weighted soil volumetric water content
#' 
#' @description Using ERA5-land SWC and a given depth, the function weight the
#' different horizons according to their thickness
#' @note ERA5-land data needs to be downloaded prior to applying the function
#' @param SWCtot dataframe of SWCtot cropped to the accurate extent
#' @param depth numeric indicating to which depth swc is to be considered
#' @return dataframe with weighted swc over different horizons
#'  
weight_swc <- function(swc,
                       depth){
  swc<-rast(swc,crs="epsg:4326")
  
  if(depth>100){
    x1=7
    x2=21
    x3=72
    x4=depth-100
  }
  if(depth<=100&depth>28){
    x1=7
    x2=21
    x3=depth-28
    x4=0
  }
  if(depth<=28&depth>7){
    x1=7
    x2=depth-7
    x3=0
    x4=0
  }
  if(depth<=7){
    x1=depth
    x2=0
    x3=0
    x4=0
  }
  swc_weighted <- (x1*swc[[grepl("h1", names(swc))]]+
                     x2*swc[[grepl("h2", names(swc))]]+
                     x3*swc[[grepl("h3", names(swc))]]+
                     x4*swc[[grepl("h4", names(swc))]])/(x1+x2+x3+x4)
  names(swc_weighted)="swc"
  return(as.data.frame(swc_weighted,xy=TRUE))
}



#' Compute psi with sureau module
#' 
#' @description function calls swc over different horizons, and use sureau model
#' to integrate their difference in conductivity in the overall psi. It take into 
#' account different possible root depth profile, maximum root depth and parameters
#' of sureau
#' @param swc 5th percentile of swc time serie, for "obs" horizons
#' @param europe Europe spatial extent
#' @param dir.hydro directory of hydraulic parameters
#' @param depth NULL if varying according to ESDAC data, numeric if set to a max
#' @param dir.depth directory of depth raster
#' @param dir.ecoregions directory of WWF ecoregions
#' @param LAImax SUREAU param: max LAI, varying from 2 to 6
#' @param fRootToLeaf SUREAU param: Root to leaf ratio
#' @param rootRadius SUREAU param: root radius 
#' @param beta SUREAU param: root profile, NULL if varying with ecoregions, or
#'  numeric if set constant 
#' @param obs vector of horizons depth of swc file
#' @param ref vector of horizons depth of hydraulic param
##' @param max_depth maximum depth considered, to extend last swc horizon
#' @param file.output file path where to write result
#' @return Psi_min dataframe, that contains values of Psi_min computed by weighting
#' the different horizons


compute_psi_sureau <- function(swc,
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
                               file.output
){
  # Rast model for resolution issues
  rast.mod<-crop(rast(paste0(dir.hydro,
                             "MRC_alp_sl2.tif")),
                 vect(europe),
                 mask=TRUE)
  europe<-vect(europe)
  
  
  # depth used
  print("Load depth")
  if (is.null(depth_max)){ # varying depth
    rast.depth=rast(dir.depth)
    crs(rast.depth) <- "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs +type=crs"
    rast.depth <- terra::project(rast.depth,"epsg:4326",method="near")
    rast.depth <- resample(rast.depth,
                           rast.mod,
                           method="near")
    names(rast.depth)="depth"
  } else{
    rast.depth=rast(rast.mod,nlyr=1,vals=depth_max)
    names(rast.depth)="depth"
  }
  
  
  ## Look for biome beta
  print("Compute betarootprofile")
  if (is.null(beta)){
    ecoregions=read_sf(dsn=dir.ecoregions, #read ecoregions
                       layer="wwf_terr_ecos") |> 
      select(BIOME,geometry)
    sf_use_s2(FALSE)
    ecoregions=st_crop(ecoregions,sf::st_bbox(europe))
    points=st_as_sf(swc[,c("x","y")],coords=c("x","y"),crs=st_crs(ecoregions))
    points_biome=st_join(points,ecoregions) |> # associate each ecoR to beta
      mutate(beta=case_when(BIOME%in%c(4,8,9,98)~0.966,
                            BIOME%in%c(5)~0.976,
                            BIOME%in%c(6)~0.943,
                            BIOME==11~0.914,
                            BIOME==12~0.964))
    beta=resample(rast(cbind(swc[,c("x","y")],
                             betaRootProfile=points_biome$beta)),
                  rast.mod,
                  method="near")
    rm(ecoregions,points,points_biome)
  } else{
    beta=resample(rast(cbind(swc[,c("x","y")],
                             betaRootProfile=beta)),
                  rast.mod,
                  method="near")
  }
  
  ## Compute psi and ksoil and ksoiltostem ##
  print("compute SUREAU module")
  # "downscale" SWC 
  swc=resample(rast(swc,crs="epsg:4326"),
               rast.mod,
               method="near")
  
  # Compute root density characteristics per horizons according to depth
  compute_horizon_density <- function(d, betaRootProfile, horizon_start, horizon_end) {
    # Use ifelse for a vectorized condition
    result <- ifelse(d <= horizon_start, 0,
                     switch(
                       as.character(horizon_end),
                       "5" = ifelse(d < 5, 1, (1 - betaRootProfile^5) / (1 - betaRootProfile^d)),
                       ifelse(d<horizon_end,
                              (betaRootProfile^horizon_start - betaRootProfile^d) / (1 - betaRootProfile^d),
                              (betaRootProfile^horizon_start - betaRootProfile^horizon_end) / (1 - betaRootProfile^d)
                       )
                     )
    )
    return(result)
  }
  
  ## Compute SUREAU ##
  RAI = LAImax*fRootToLeaf
  compute.B_GC <- function(hor_value,
                           depth,
                           RAI, 
                           rootRadius,
                           horizon_start,
                           horizon_end
  ) {
    d=ifelse(depth<=horizon_start,
             NA,
             ifelse(depth>horizon_end,
                    horizon_end,
                    depth))
    La <- RAI * hor_value / (2 * pi * rootRadius)
    Lv <- RAI * hor_value / (2 * pi * rootRadius * d)
    b <- 1 / sqrt(pi * Lv)
    B_GC <- La * 2 * pi / log(b / rootRadius)
    return(B_GC)
  }
  
  rast.depth.samp <- c(rast.depth, beta) %>%
    as.data.frame(xy = TRUE, na.rm = FALSE) %>%
    mutate(depth = as.numeric(na_if(depth, 0))) %>%
    mutate(
      hor.1 = compute_horizon_density(d=depth, betaRootProfile, 0, 5),
      hor.2 = compute_horizon_density(d=depth, betaRootProfile, 5, 15),
      hor.3 = compute_horizon_density(d=depth, betaRootProfile, 15, 30),
      hor.4 = compute_horizon_density(d=depth, betaRootProfile, 30, 60),
      hor.5 = compute_horizon_density(d=depth, betaRootProfile, 60, 100),
      hor.6 = compute_horizon_density(d=depth, betaRootProfile, 100, 200),
      hor.7 = compute_horizon_density(d=depth, betaRootProfile, 200, Inf)
    ) |> 
    # compute bcg per horizon
    mutate(hor.1=compute.B_GC(hor.1,depth,RAI,rootRadius,0,5),
           hor.2=compute.B_GC(hor.2,depth,RAI,rootRadius,5,15),
           hor.3=compute.B_GC(hor.3,depth,RAI,rootRadius,15,30),
           hor.4=compute.B_GC(hor.4,depth,RAI,rootRadius,30,60),
           hor.5=compute.B_GC(hor.5,depth,RAI,rootRadius,60,100),
           hor.6=compute.B_GC(hor.6,depth,RAI,rootRadius,100,200),
           hor.7=compute.B_GC(hor.7,depth,RAI,rootRadius,200,Inf)) |> 
    select(-betaRootProfile)
  rast.depth.samp=rast(rast.depth.samp,crs="epsg:4326")
  rm(rast.depth,beta)
  #rast(SWCmin,crs="epsg:4326")
  
  
  # correspondance between horizons fril ERA5 and 3D hydrosoil
  calculate_weights <- function(ref, obs,max_depth=3) {
    
    if(ref[length(ref)]<max_depth){
      ref_ext  <- c(ref,max_depth)
    } else(ref_ext  <- ref)
    if(obs[length(obs)]<ref_ext[length(ref_ext)]){
      obs_ext  <- c(obs,ref_ext[length(ref_ext)])
    } else(obs_ext  <- obs)    
    
    # Create a matrix filled with zeros
    weight_matrix <- matrix(0, nrow =length(obs_ext)-1, ncol = length(ref_ext)-1)
    
    for(i in 1:(length(obs_ext)-1)) {
      for(j in 1:(length(ref_ext)-1)) {
        if(obs_ext[i]<ref_ext[j+1] & ref_ext[j]<obs_ext[i+1]){
          weight_matrix[i,j]=min(1,
                                 (obs_ext[i+1]-ref_ext[j])/(ref_ext[j+1]-ref_ext[j]),
                                 (ref_ext[j+1]-obs_ext[i])/(ref_ext[j+1]-ref_ext[j]),
                                 (obs_ext[i+1]-obs_ext[i])/(ref_ext[j+1]-ref_ext[j]),
                                 na.rm=TRUE)
        }
      }
    }
    if (sum(colSums(weight_matrix))!=length(ref)){print("error")}
    return(weight_matrix)
  }
  
  weight_matrix<-calculate_weights(ref,
                                   obs)
  
  # add a fictive layer if deepest horizon doesn't reach max_depth
  if(obs[length(obs)]<max(ref[length(ref)],max_depth)){
    swc=c(swc,swc[[nlyr(swc)]])
  }
  
  swc_weighted=rast(swc,nlyr=length(ref),vals=NA)
  
  
  print("Computing weighted swc")
  for(i in 1:nlyr(swc_weighted)){
    print(i)
    swc_weighted[[i]]=sum(weight_matrix[,i]*swc,na.rm = TRUE)
  }
  
  
  
  # Gather spatialized PDF parameters with SWC min/max for each of the 7 horizons
  
  print("Croping hydraulic pars")
  pars.files=list.files(dir.hydro)
  mrc.files=pars.files[grepl("MRC_",pars.files)]
  mrc.rast=crop(rast(file.path(dir.hydro,mrc.files)),
                europe)
  hcc.files=pars.files[grepl("HCC_",pars.files)]
  hcc.rast=crop(rast(file.path(dir.hydro,hcc.files)),
                europe)  
  names(swc_weighted)=rep("swc_min",length(ref))
  names(rast.depth.samp)=c("depth",rep("BGC",length(ref)))
  
  
  print("Compute psi per horizon")
  psi_hor=as.data.frame(swc_weighted,xy=TRUE,na.rm=FALSE) |> select(x,y)
  for (i in 1:7){ # loop on the 7 horizons of terraclimate
    print(i)
    psi=c(swc_weighted[[i]],
          rast.depth.samp[[1]], # depth 
          rast.depth.samp[[i+1]], # BCG of horizons
          mrc.rast[[grep(paste0("sl",i),names(mrc.rast))]],
          hcc.rast[[grep(paste0("sl",i),names(hcc.rast))]]
    )
    names(psi)=str_replace(names(psi), paste0("_sl",i), "") #make generic names
    
    psi=as.data.frame(psi,xy=TRUE,na.rm=FALSE)
    setDT(psi)
    psi[, ':=' (
      REW_mrc = ifelse((swc_min - MRC_thr*10^(-4)) / ((MRC_ths - MRC_thr) * 10^(-4)) < 0, 0.01,
                       ifelse((swc_min - MRC_thr*10^(-4)) / ((MRC_ths - MRC_thr) * 10^(-4)) > 1, 1,
                              (swc_min - MRC_thr*10^(-4)) / ((MRC_ths - MRC_thr) * 10^(-4)))),
      REW_hcc = ifelse((swc_min - HCC_thr*10^(-4)) / ((HCC_ths - HCC_thr) * 10^(-4)) < 0, 0.01,
                       ifelse((swc_min - HCC_thr*10^(-4)) / ((HCC_ths - HCC_thr) * 10^(-4)) > 1, 1,
                              (swc_min - HCC_thr*10^(-4)) / ((HCC_ths - HCC_thr) * 10^(-4)))))]
    psi[, ':=' (
      psi_min = -(((1/REW_mrc)^(1/(MRC_m*10^(-4))) - 1)^(1/(MRC_n*10^(-4))) * (1/(MRC_alp*10^(-4))) * 9.78*10^(-2)),
      ksoil = HCC_K0 * (REW_hcc^(HCC_L*10^(-4))) * (1 - (1 - REW_hcc^(1/(HCC_m*10^(-4))))^(HCC_m*10^(-4)))^2,
      ksoilGC = (HCC_K0 * (REW_hcc^(HCC_L*10^(-4))) * (1 - (1 - REW_hcc^(1/(HCC_m*10^(-4))))^(HCC_m*10^(-4)))^2) * 1000 * BGC)]
    psi[, ':=' (
      psi_w = ksoilGC * psi_min
    )]
    colnames(psi)[!colnames(psi) %in% c("x","y","depth")]=paste0(colnames(psi)[!colnames(psi) %in% c("x","y","depth")],"_",i)
    psi_hor=cbind(psi_hor,
                  swc_min=psi$swc_min,
                  ksoilGC=psi$ksoilGC,
                  psi_w=psi$psi_w) |> 
      rename(!!paste0("swc_min_",i):="swc_min",
             !!paste0("ksoilGC_",i):="ksoilGC",
             !!paste0("psi_w_",i):="psi_w") 
    rm(psi)
  }
  rm(hcc.rast,mrc.rast)
  
  print("Compute psi total")
  all_psi_summary = psi_hor %>%
    select(x, y, matches("ksoilGC_|psi_w_")) %>%
    filter(!is.na(psi_w_1)) 
  
  all_psi_summary$sum_psi=rowSums(all_psi_summary[grepl("psi_w_",colnames(all_psi_summary))],na.rm = TRUE)
  all_psi_summary$sum_k=rowSums(all_psi_summary[grepl("ksoilGC_",colnames(all_psi_summary))],na.rm = TRUE)
  all_psi_summary$psi=all_psi_summary$sum_psi/all_psi_summary$sum_k 
  
  fwrite(all_psi_summary,file.output)
  
  return(all_psi_summary)
}


#'Get LT50 database
#'
#'@description From Constance database and measures performed in Clermont campaign
#'2022, we build a dataset with consistent measures
#'@return Dataframe with LT50 per European species and sd associated when available

get_LT50 <- function(file.output){
  ### Load LT50 database from litterature ###
  ###########################################
  df.LT50.lit <- read.csv2(file = "data/Species traits/df_LT50_sp_checked.csv", header = T,na.strings = c("", "NA")) # made by code in "R/matchsp_get_tax.R"
  df.LT50.lit$New_continent <- "Continent"
  # remove NA lines (where temperature was reported as "below xx°")
  df.LT50.lit %>% filter(!is.na(Temperature_of_resistance)) -> df.LT50.lit 
  
  ### cleaning data coding ###
  ############################
  df.LT50.lit %>%
    mutate(across(.cols=c("Freezing_rate_.K.h.1.","Temperature_of_resistance","Duration_of_min_temp_exposure_.hours.","Longitude","Latitude"),
                  ~as.numeric(.))) %>% 
    mutate(Country = case_when(!is.na(Country) ~ Country,
                               # Morin et al Tree Physiol 2007
                               Reference == "Morin et al Tree Physiol 2007" & Species == "Quercus robur" ~ "Germany",
                               Reference == "Morin et al Tree Physiol 2007" & Species == "Quercus ilex" ~ "France",
                               Reference == "Morin et al Tree Physiol 2007" & Species == "Quercus pubescens" ~ "France",
                               # Sakai et al Ecology 1981 correct country below
                               Reference == "Sakai et al Ecology 1981" & Location %in% c("Canberra Botanical Garden", "Craigieburn", "Mount Ginini", "Snowy Mountains", "Tasmania") ~ "Australia",
                               Reference == "Sakai et al Ecology 1981" & Location %in% c("Christchurch Botanical Garden", "Arthurs Pass", "Jacksons", "Porters Pass", "Otira", "Waimakariri") ~ "New Zealand",
                               Reference == "Sakai et al Ecology 1981" & Location %in% c("Mount Fuji", "Hakodate", "Mount Tsukuba", "Sapporo, Hokkaido", "Shiojiri", "Tokyo", "Yamabe, Hokkaido") ~ "Japan",
                               Reference == "Sakai et al Ecology 1981" & Species %in% c("Araucaria heterophylla", "Callitris endlicheri", "Diselma archeri",  "Microcachrys tetragona") ~ "Australia", 
                               Reference == "Sakai et al Ecology 1981" & Species %in% c("Callitris oblonga", "Dacrycarpus dacrydioides", "Dacrydium cupressinum", "Halocarpus bidwillii", "Libocedrus bidwillii", "Phyllocladus alpinus", "Pinus pseudostrobus", "Podocarpus nivalis") ~ "New Zealand",
                               # Bannister New Zealand J. of Bot. 1990 al grown in NZ
                               Reference == "Bannister New Zealand J. of Bot. 1990" ~ "New Zealand",
                               # Sakai Can J. Bot. 1983 "himalaya" species were sent from Bot garden in PNG or from field in Nepal
                               Reference == "Sakai Can J. Bot. 1983" & Location == "Himalaya" & Species %in% c("Abies spectabilis", "Larix potaninii", "Tsuga dumosa", "Juniperus monosperma", "Juniperus squamata", "Picea smithiana") ~ "Nepal",
                               Reference == "Sakai Can J. Bot. 1983" & Species %in% c("Picea orientalis") ~ "Japan",
                               Reference == "Sakai Can J. Bot. 1983" & Species %in% c("Picea smithiana") ~ "Nepal",
                               Reference == "Darrow et al New Zealand J. of Bot. 2001" ~ "New Zealand",
                               Comment == "Reported from Oohata & Sakai 1982" ~ "Japan",
                               Reference == "Read & Hill Aust. J. Bot. 1988" ~ "Australia",
                               Reference == "Harrison et al Plant Physiol. 1978" ~ "USA",
                               Reference == "Alberdi et al Phytochemistry 1989" ~ "Chile",
                               Reference ==  "Bannister New Zealand J. of Bot. 2003" ~ "New Zealand",
                               Reference ==  "Goldstein et al Oecologia 1985" ~ "Venezuela",
                               TRUE ~ Country),
           Organ = case_when(!is.na(Organ) ~ Organ,
                             Comment == "Reported in Bannister New Zealand J. of Bot. 2007" ~ "Leaf",
                             Reference == "Read & Hill Aust. J. Bot. 1988" ~ "Leaf",
                             Reference == "Durham et al Physiol. Plant. 1991" ~ "Leaf"),
           Freezing_rate_.K.h.1. = case_when(!is.na(Freezing_rate_.K.h.1.) ~ Freezing_rate_.K.h.1.,
                                             Reference == "Darrow et al New Zealand J. of Bot. 2001" ~ 4,
                                             Reference ==  "Goldstein et al Oecologia 1985" ~ 10,
                                             TRUE ~ Freezing_rate_.K.h.1.),
           Method = case_when(!is.na(Method) ~ Method,
                              Reference == "Darrow et al New Zealand J. of Bot. 2001" ~ "Visual scoring",
                              Reference == "Alberdi et al Phytochemistry 1989" ~ "Visual scoring",
                              Reference ==  "Bannister New Zealand J. of Bot. 2003" ~ "Visual scoring",
                              Reference ==  "Goldstein et al Oecologia 1985" ~ "TTC",
                              TRUE ~ Method)
    ) -> df.LT50.lit
  
  ### Load LT50 database from 2022 campaign ###
  #############################################
  df.LT50.2022 <- read.csv2("data/Species traits/camp_LT50_DB.csv") |> # made by code in "Analysis new campaign data.Rmd" : output file
    bind_rows(read.csv2("data/Species traits/new_data_2023_prelim.csv"))
  
  
  #clean rm species with bad fits:
  df.LT50.2022 %>% filter(!Species %in% c("Abies balsamea", 
                                          "Acer saccharum",
                                          "Betula papyrifera",
                                          "Larix laricina",
                                          "Tsuga canadensis", 
                                          "Pinus banksiana",
                                          "Picea mariana")) -> df.LT50.2022
  df.LT50.2022$Thawing_rate_.K.h.1. <- as.character(df.LT50.2022$Thawing_rate_.K.h.1.)
  
  ### merge together ###
  ######################
  full_join(df.LT50.lit, df.LT50.2022) -> df.LT50
  df.LT50 <- df.LT50 %>% select(-Family, -Genus, -Phyllum)
  
  # country lists
  S_A <- c("Argentina", "Brazil", "Chile", "South America", "Venezuela")
  Asia <- c("China", "Iran", "Israel" , "Japan", "Korea", "Russia", "Taiwan", "Turkey", "Asia Minor", "Himalaya", "Nepal")
  E_U <- c("Austria", "Croatia", "Czech Rep.", "Denmark", "England", "Finland", "France", "Germany", "Iceland", "Italy", "Poland", "Romania", "Serbia", "Slovakia", "Spain", "Sweden", "Swiss", "Switzerland", "Ukraine", "Europe")
  N_A <- c("Canada", "Mexico", "Quebec", "USA")
  Oce <- c("Australia", "New Guinea", "Papua New Guinea", "New Zealand")
  Africa <- c("South Africa", "Canary Islands")
  
  df.LT50 %>% mutate(New_continent = case_when(Country %in% E_U  | Provenance %in% E_U ~ "E_U",
                                               Country %in% Asia | Provenance %in% Asia  ~ "Asia",
                                               Country %in% N_A | Provenance %in% N_A  ~ "N_A",
                                               Country %in% S_A | Provenance %in% S_A   ~ "S_A",
                                               Country %in% Africa | Provenance %in% Africa  ~ "Africa",
                                               Country %in% Oce | Provenance %in% Oce ~ "Oceania",
                                               is.na(Country) ~ Provenance,
                                               TRUE ~ Provenance)) -> df.LT50
  rm(S_A,Africa,Asia,E_U,N_A,Oce)
  
  bud <- c("Lateral bud",  "Bud", "Bud base vascular tissue", "Leaf primordia", "Primordial shoot", "Floral primordial", "Procambium")
  flowerbud <- c("Flower Bud", "Flower bud", "Female flower", "Male flower")
  branch <- c("Cane", "Basal stem","Wood", "Twig", "Shoot", "Stem", "Upper-crown shoot", "Twig cambium", "Xylem", "Xylem parenchyma", "Shoot vascular tissue", "Pith", "Pith parenchyma", "Phloem", "Cambial meristem", "Cambium", "Cortex", "Branch", "branch")
  leaf <- c("Needle", "Inner crown leaf", "Leaf", "Upper crown leaf")
  root <- c("Root", "Root cambium")
  df.LT50 <- df.LT50 %>% filter(Organ %in% c(bud, flowerbud, branch, leaf, root)) 
  
  # classify in organ type
  df.LT50$OrganGroup <- NA
  df.LT50$OrganGroup[df.LT50$Organ %in% bud] <- "bud"
  df.LT50$OrganGroup[df.LT50$Organ %in% flowerbud] <- "flowerbud"
  df.LT50$OrganGroup[df.LT50$Organ %in% branch] <- "branch"
  df.LT50$OrganGroup[df.LT50$Organ %in% leaf] <- "leaf"
  df.LT50$OrganGroup[df.LT50$Organ %in% root] <- "root"
  
  rm(bud,flowerbud,branch,leaf,root)
  
  names(df.LT50) <- gsub(names(df.LT50), pattern = " ", replacement = "_")
  
  ### filter for analysis ###
  ###########################
  
  growthform.select <- c("tree","Tree") #"Small tree",
  plantage.select <- c("Mature",NA)
  organgroup.select <-c("branch","bud","flowerbud","leaf") #
  method.select <- c("EL", "Visual scoring") #
  period <- c("Winter")
  
  df.LT50 %>% 
    filter(Growth_form %in% growthform.select &
             Period %in% period &
             #Plant_age %in% plantage.select &
             OrganGroup %in% organgroup.select &
             Method %in% method.select 
    ) %>%
    filter(Temperature_of_resistance>(-100)) %>% 
    filter(!(Species %in% c("Abies balsamea",
                            "Betula platyphylla",
                            "Salix sachalinensis", 
                            "Larix laricina", 
                            "Pinus banksiana", 
                            "Pinus strobus", 
                            "Thuja occidentalis",
                            "Pinus resinosa") &
               Temperature_of_resistance == -196)) %>%
    filter(!((Species=="Larix decidua")& (Reference=="Charrier et al Tree Phys 2013"))) %>% # measure that do not exist
    filter(!(Reference=="Vitra et al New Physiol. 2017" & Date!="day 29")) %>% #remove measures coming from temporal dehardening series
    filter(!(Species %in% c("Juglans nigra") &
               Organ == "Cortex")) %>% 
    filter(Type == "LT50") %>% 
    filter(New_continent=="E_U") %>% 
    mutate(data.quality=NA) -> df.LT50.select
  
  ### Test of filtering ###
  #########################
  unique(df.LT50.select$Thawing_rate_.K.h.1.)
  df.LT50.select %>% 
    filter(Plant_age == "Mature" & 
             Fit %in% c("sigmoid","Sigmoid") &
             Method == "EL" &
             OrganGroup == "branch") %>%     
    group_by(Reference,Species) %>%
    slice(which.min(Temperature_of_resistance)) %>%
    ungroup()  -> df.LT50.select.strong
  
  
  thawing.rate <- c("2","3","5","7")
  #sp.deleted <- setdiff(unique(df.LT50.select$Species),unique(df.LT50.select.strong$Species))
  df.filtered <- as.data.frame(matrix(nrow = 0,
                                      ncol=dim(df.LT50.select)[2]))
  colnames(df.filtered) <- colnames(df.LT50.select)
  for (sp in unique(df.LT50.select$Species)){
    print(sp)
    print("step 1")
    #### All constraints
    df.sp <- df.LT50.select[df.LT50.select$Species==sp & df.LT50.select$Plant_age == "Mature" & 
                              df.LT50.select$Thawing_rate_.K.h.1. %in% thawing.rate &
                              df.LT50.select$Method == "EL" & 
                              df.LT50.select$OrganGroup == "branch" &
                              df.LT50.select$Fit %in% c("sigmoid","Sigmoid")
                            ,] %>%     
      group_by(Reference,Species) %>%
      slice(which.min(Temperature_of_resistance)) %>%
      ungroup() 
    if (dim(df.sp)[1]>0){
      df.sp$data.quality=1
      df.filtered <- rbind(df.filtered,df.sp)
    }else{
      #### Release "sigmoid" constraint
      print("step 2")
      df.sp <- df.LT50.select[df.LT50.select$Species==sp & df.LT50.select$Plant_age == "Mature" & 
                                df.LT50.select$Thawing_rate_.K.h.1. %in% thawing.rate &
                                df.LT50.select$Method == "EL" & 
                                df.LT50.select$OrganGroup == "branch" 
                              ,] %>%     
        group_by(Reference,Species) %>%
        slice(which.min(Temperature_of_resistance)) %>%
        ungroup() 
      if (dim(df.sp)[1]>0){
        df.sp$data.quality=2
        df.filtered <- rbind(df.filtered,df.sp)
      }else{
        #### Release "branch" only constraint t "branch" "bud"
        print("step 3")
        df.sp <- df.LT50.select[df.LT50.select$Species==sp & df.LT50.select$Plant_age == "Mature" & 
                                  df.LT50.select$Thawing_rate_.K.h.1. %in% thawing.rate &
                                  df.LT50.select$Method == "EL" & 
                                  df.LT50.select$OrganGroup %in% c("branch","bud")
                                ,] %>%     
          group_by(Reference,Species) %>%
          slice(which.min(Temperature_of_resistance)) %>%
          ungroup() 
        if (dim(df.sp)[1]>0){
          df.sp$data.quality=3
          df.filtered <- rbind(df.filtered,df.sp)
        }else{
          #### Release "method" constraint
          print("step 4")
          df.sp <- df.LT50.select[df.LT50.select$Species==sp & df.LT50.select$Plant_age == "Mature" &
                                    df.LT50.select$Thawing_rate_.K.h.1. %in% thawing.rate &
                                    df.LT50.select$OrganGroup %in% c("branch","bud")
                                  ,]%>%     
            group_by(Reference,Species) %>%
            slice(which.min(Temperature_of_resistance)) %>%
            ungroup() 
          if (dim(df.sp)[1]>0){
            df.sp$data.quality=4
            df.filtered <- rbind(df.filtered,df.sp)
          }else{
            #### Release "thawing rate" constraint
            print("step 5")
            df.sp <- df.LT50.select[df.LT50.select$Species==sp & df.LT50.select$Plant_age == "Mature" &
                                      df.LT50.select$OrganGroup %in% c("branch","bud")
                                    ,]%>%     
              group_by(Reference,Species) %>%
              slice(which.min(Temperature_of_resistance)) %>%
              ungroup() 
            if (dim(df.sp)[1]>0){
              df.sp$data.quality=5
              df.filtered <- rbind(df.filtered,df.sp)
            }else{
              #### Release "age" constraint
              print("step 6")
              df.sp <- df.LT50.select[df.LT50.select$Species==sp & df.LT50.select$OrganGroup %in% c("branch","bud")
                                      ,] %>%     
                group_by(Reference,Species) %>%
                slice(which.min(Temperature_of_resistance)) %>%
                ungroup() 
              if (dim(df.sp)[1]>0){
                df.sp$data.quality=6
                df.filtered <- rbind(df.filtered,df.sp) 
              }else{
                print("step 7")
                df.filtered <- rbind(df.filtered,
                                     df.LT50.select[df.LT50.select$Species==sp
                                                    ,]%>%     
                                       group_by(Reference,Species) %>%
                                       slice(which.min(Temperature_of_resistance)) %>%
                                       ungroup() %>% 
                                       mutate(data.quality=7))
              }
            }
          }
        }
      }
    }
    rm(df.sp)}
  
  
  ### Summarise per species ###
  #############################
  df.species.traits <- df.filtered %>% 
    filter(New_continent=="E_U") %>%
    group_by(Species,Period,data.quality) %>% 
    summarise(LT50.mean=mean(Temperature_of_resistance,na.rm=TRUE),
              LT50.sd=sd(Temperature_of_resistance,na.rm=TRUE)) %>% 
    ungroup() %>% 
    arrange(LT50.mean) %>% 
    mutate(rank_filt=row_number())
  
  ## For raw data
  df.species.traits.raw <- df.LT50.select %>% 
    filter(New_continent=="E_U") %>%
    group_by(Species,Period) %>% 
    summarise(LT50.mean=mean(Temperature_of_resistance,na.rm=TRUE),
              LT50.sd=sd(Temperature_of_resistance,na.rm=TRUE)) %>% 
    ungroup() %>% 
    arrange(LT50.mean) %>% 
    mutate(rank=row_number()) %>% 
    left_join(df.species.traits[,c("Species","rank_filt")],by="Species")
  
  
  
  
  ### Few plots ###
  #################
  # ## Plot difference between mean computed over raw dataset and with filtered dataset
  # df.species.traits %>%
  #   filter(!is.na(LT50.sd)) %>%
  #   rename(`LT50.mean.filter`="LT50.mean",
  #          `LT50.sd.filter`="LT50.sd") %>%
  #   left_join(df.species.traits.raw %>%
  #               filter(!is.na(LT50.sd)) %>%
  #               rename(`LT50.mean.raw`="LT50.mean",
  #                      `LT50.sd.raw`="LT50.sd"),
  #             by="Species") %>%
  #   ggplot(aes(LT50.mean.filter-LT50.mean.raw,Species))+
  #   geom_point()
  # 
  # ## Boxplots for species with several measures
  # list=df.filtered  %>%
  #   mutate(data.quality=as.factor(data.quality)) %>% 
  #   filter(New_continent=="E_U") %>%
  #   # group_by(Species) %>%
  #   # filter(n()>1) %>%
  #   # ungroup() %>%
  #   ggplot(aes(Temperature_of_resistance,Species))+
  #   geom_line()+
  #   geom_point(aes(shape=data.quality))+
  #   theme_bw()+
  #   theme(axis.title=element_blank() )
  # ggsave(list,
  #        filename = "species_list.pdf",
  #        path="output/",
  #        device="pdf",
  #        height=12)
  # 
  # ## Plots with changes in species ranking before/after LT50 filtering
  # df.species.traits.raw %>% 
  #   filter(abs(rank-rank_filt)>5) %>% 
  #   ggplot(aes(rank-rank_filt,Species))+
  #   geom_point()
  # 
  # 
  # df.species.traits.raw %>%
  #   pivot_longer(cols = c("rank","rank_filt"),names_to = "rank.type",values_to = "rank") %>%
  #   mutate(rank.type=as.factor(rank.type)) %>%
  #   ggplot(aes(x = rank.type, y = rank, group = Species)) +
  #   geom_line(aes(color = Species, alpha = 1), size = 2) +
  #   geom_point(aes(color = Species, alpha = 1), size = 2) +
  #   geom_point(color = "#FFFFFF", size = 1) +
  #   scale_y_reverse(breaks = 1:nrow(df.species.traits.raw)) +
  #   scale_x_discrete(breaks = 1:10) +
  #   theme_bw()+
  #   theme(legend.position = 'none',
  #         axis.title = element_blank(),
  #         axis.text.y = element_blank(),
  #         axis.ticks.y= element_blank()) +
  #   geom_text(data = df.species.traits.raw %>% arrange(rank),
  #             aes(label = Species, x = .95) , hjust = .5,
  #             color = "#888888", size = 3) +
  #   labs(x = '', y = 'Rank', title = 'Changes in LT50 species ranking after filtering')
  # 
  # 
  
  
  write.csv(df.species.traits,file.output,row.names = FALSE)
  
  return(list(df.LT50.raw=df.LT50.select,
              df.LT50.cor=df.filtered,
              df.LT50sp.raw=df.species.traits.raw,
              df.LT50sp.cor=df.species.traits))
  
}




#' Get P50 database 
#' 
#' @description Load P50 from Max database
#' @return  Dataframe with P50 per European species 
get_P50 <- function(species.list,
                    file.output){
  df.p50.msp<-read.csv2("data/Species traits/p50_nmsp.csv") |> 
    filter(species.binomial %in% species.list) |> 
    mutate(group=tolower(group),
           p88.mu=p50-50/slope,
           p50.mu=p50,
           bdd="martinsp") |> 
    select(group,species.binomial,p50.mu,p88.mu,slope,bdd) |> 
    filter(!(group=="angiosperm"&is.na(p88.mu)))
  df.p50.lit <- read.csv2("data/Species traits/2022_10_20_cleaned_xft.csv") %>% 
    select(uniqueID,lat,long,XFT.database,Group,Family,Genus,Species,
           Developmental.stage,Growth.form,P50,P50.SD,P12,P88,Curve,Equation,
           psi.min.predawn..MPa.,psi..min.midday..MPa.,P50.number.of.samples) %>% 
    mutate(species.binomial=paste(Genus,Species),
           bdd="hammond") %>% 
    filter((species.binomial %in% species.list)&
             (!species.binomial %in% df.p50.msp$species.binomial )) |> 
    filter(Curve=="S") %>%
    filter(Growth.form=="T") %>% 
    filter(Developmental.stage =="A") |> 
    group_by(Group,Genus,Species,species.binomial,bdd) %>% 
    summarise(P50.mu=mean(as.numeric(P50),na.rm=TRUE),
              P50.sd=sd(as.numeric(P50),na.rm=TRUE),
              P88.mu=mean(as.numeric(P88),na.rm=TRUE),
              P88.sd=sd(as.numeric(P88),na.rm=TRUE)) %>% 
    ungroup() |> 
    rename_with(.cols=everything(),
                tolower) |> 
    select(group,species.binomial,p50.mu,p50.sd,p88.mu,p88.sd,bdd)
  df.p50=bind_rows(df.p50.lit,df.p50.msp)
  write.csv(df.p50,file.output,row.names = FALSE)
  return(df.p50)
}

#' Get LT50/P50 database 
#' 
#' @description Load LT50/P50 for european species
#' @param dir.data data directory
#' @param dir.file directory of species file
#' @return dataframe of LT50/P50 traits per species
#' 

get_traits <- function(dir.distribution="data/chorological_maps_dataset",
                       df.P50, #=df.P50
                       df.LT50,
                       file.output){ #=df.LT50$df.LT50sp.cor
  # df.P50 <- read.csv("output/df_P50_filtered.csv") 
  # df.LT50 <- read.csv("output/df_LT50_filtered.csv")
  
  df.traits.raw <- df.LT50 %>% 
    left_join(df.P50,by=c("Species"="species.binomial")) |> 
    rename_with(.cols=everything(),
                tolower)
  # df.traits.raw <- full_join(df.P50,df.LT50,by="Species")
  # df.traits.missing <- df.traits.raw %>% 
  #   filter(is.na(P50tot) | is.na(LT50.mean))
  df.traits <- df.traits.raw %>% 
    filter(!is.na(p50.mu)) %>% 
    filter(!is.na(lt50.mean)) %>% 
    select(species,group,p50.mu,p50.sd,p88.mu,p88.sd,slope,lt50.mean,lt50.sd,data.quality) %>% #MATmean,MAPmean,Leaf_phenology,
    unique() %>% 
    rename(`species.binomial`="species") %>% 
    separate(col = species.binomial,
             into=c("genus","species"),
             remove=FALSE) %>% 
    mutate(sp.ind=paste0(substr(genus,1,2),substr(species,1,2)),
           species.name=str_replace(species.binomial," ","_")) %>% 
    relocate(c("species.name","sp.ind"),.after="species") %>% 
    mutate(p.trait=case_when(group=="gymnosperm"~"p50",
                             group=="angiosperm"&!is.na(p88.mu)~"p88",
                             group=="angiosperm"&is.na(p88.mu)~"p50"),
           px.mu=case_when(p.trait=="p88"~p88.mu,
                           p.trait=="p50"~p50.mu),
           px.sd=case_when(p.trait=="p88"~p88.sd,
                           p.trait=="p50"~p50.sd)) %>% 
    filter(!(species.binomial=="Pinus contorta"&group=="angiosperm"))
  
  for (i in 1:dim(df.traits)[1]){
    species.files=list.files(file.path(dir.distribution,df.traits$species.binomial[i],"shapefiles"))
    if(length(species.files[grepl(paste0(df.traits$species.name[i],"_plg_clip"),species.files)])>1){
      df.traits$file[i]=paste0(df.traits$species.name[i],"_plg_clip")
    } else if(length(species.files[grepl(paste0(df.traits$species.name[i],"_",df.traits$species[i],"_plg_clip"),species.files)])>1){
      df.traits$file[i]=paste0(df.traits$species.name[i],"_",df.traits$species[i],"_plg_clip")
    } else if(length(species.files[grepl(paste0(df.traits$species.name[i],"_plg"),species.files)])>1) {
      df.traits$file[i]=paste0(df.traits$species.name[i],"_plg")
    } else if(length(species.files[grepl(paste0(df.traits$species.name[i],"_",df.traits$species[i],"_plg"),species.files)])>1) {
      df.traits$file[i]=paste0(df.traits$species.name[i],"_",df.traits$species[i],"_plg")
    } else {
      df.traits$file[i]=NA
    }
  }  
  write.csv(df.traits,file.output,row.names = FALSE)
  return(df.traits)
}



#' Get LT50/P50 database, from MAX
#' 
#' @description Load LT50/P50 for european species
#' @param dir.file directory of species file
#' @return dataframe of LT50/P50 traits per species
#' 
get_traits_max <-function(species.list,
                          dir.file="data/Species traits/base_traits_P50_LTx.xlsx",
                          file.output){
  p50.msp<-read.csv2("data/Species traits/p50_nmsp.csv") |> 
    rename(slope_msp=slope,species=species.binomial)
  traits<-readxl::read_xlsx("data/Species traits/base_traits_P50_LTx.xlsx")|>
    filter(Species%in% species.list) |> 
    rename_with(.cols=everything(),
                tolower) |> 
    left_join(p50.msp[,c("species","slope_msp")]) |>  
    mutate(taxa=case_when(class=="Magnoliopsida"~"angiosperm",
                          class=="Pinopsida"~"gymnosperm"),
           px=case_when(taxa=="angiosperm"~ifelse(!is.na(p88),p88,p50-50/slope_msp),
                        taxa=="gymnosperm"~p50),
           ptrait=case_when(taxa=="angiosperm"~"p88",
                            taxa=="gymnosperm"~"p50")
    ) |> 
    # relocate(taxa,.after=genus) |> 
    # relocate(px,.before = p50) |> 
    filter(!is.na(px)&!is.na(ltx)) |> 
    mutate(lt50_qual=case_when(is.na(ltx_clean)~"low",
                               TRUE~"ok"),
           ltx_clean_2=case_when(is.na(ltx_clean)~ltx,
                                 TRUE~ltx_clean),
           source_ltx_2=case_when(is.na(ltx_clean)~source_ltx,
                                  TRUE~source_ltx_clean)) |> 
    select(species,class,order,family,genus,taxa,
           ltx_clean_2,ltx_clean_sd,source_ltx_2,lt50_qual,
           px,p50,p50sd,p12,p88,slope,slope_msp,source) |> 
    filter(!grepl("LT0_",source_ltx_2)) |> 
    rename(
      lt50=ltx_clean_2,
      lt50_source=source_ltx_2,
      lt50_sd=ltx_clean_sd,
      px_source=source
    )
  write.csv(traits,file.output)
  return(traits)
}

