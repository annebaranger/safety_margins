#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#### SCRIPT INTRODUCTION ####
#
#' @name functions_analyses.R  
#' @description R script containing all functions relative to data
#               analyses
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Load Mauri pres/abs and format dataset
#' 
#' @description load Mauri dataset, format presence/absence and load safety margins
#' and others climatics variables associated to pres/abs points
#' @note Chelsea dwonloaded from website
#' @param dir.occ directory of mauri abs.pres dataset
#' @param df.traits dataframe of traits associated with eahc species
#' @param psi_min psimin computed over europe
#' @param frost.index frost index computed over europe
#' @return dataframe, where for each location, abs/pres of all selected species
#' is mentionned, and associated predictors

get_mauri <- function(dir.occ,species.list){
  df.Euforest=read.table(dir.occ,
                         header=TRUE,
                         sep=",",
                         dec=".")%>% 
    #remove duplicated location_tree 
    filter(!duplicated(paste(X,Y, SPECIES.NAME))) %>%
    #create plot id
    mutate(Pres = rowSums(cbind(NFI==1,FF==1,BS==1)),
           Node = unclass(factor(paste(X, Y))))
  
  db.cont=df.Euforest %>% 
    group_by(Node,SPECIES.NAME) %>%
    summarise(n=n())%>%
    spread(SPECIES.NAME, n) %>% 
    ungroup() %>%
    pivot_longer(colnames(.)[-1],
                 names_to = "species",
                 values_to = "presence") %>% 
    filter(species %in% species.list) %>% 
    left_join(df.Euforest[!duplicated(df.Euforest$Node), c("X","Y","COUNTRY","EEO","Node")],
              by = "Node") 
  
  species.list.eff = db.cont |> 
    filter(presence==1) |> 
    group_by(species) |> 
    summarise(n=n()) |> 
    filter(n>400)
  
  db.cont<-db.cont |> 
    filter(species %in% species.list.eff$species) |> 
    rename_with(.cols=everything(),
                tolower)
  
  #Change coordinates
  coordinates(db.cont) <- c("x",  "y")
  proj4string(db.cont) <- CRS("+init=epsg:3035")
  db.cont <- as.data.frame(spTransform(db.cont, CRS("+proj=longlat +datum=WGS84")))
  
  return(db.cont)
}


get.occ.clim<-function(db.mauri,
                       clim.list,
                       file.path){
  ## Load safe ty margin indicators
  clim.list=lapply(clim.list,rast)
  df.clim=lapply(seq_along(clim.list),function(x) {
    df=extract(clim.list[[x]],data.frame(x=db.mauri$x,y=db.mauri$y))[-1]
    colnames(df)=names(clim.list)[x]
    return(df)}
  )
  
  db.mauri<-cbind(db.mauri,
                  as.data.frame(df.clim))
  
  db.mauri$presence[is.na(db.mauri$presence)] <- 0
  
  fwrite(db.mauri,file.path)
  # write.csv(df.fdg.sp,"output/LT50_spring.csv")
  
  return(db.mauri)
}


#' Compute safety margins on abs/pres dataset
#' 
#' @description merge with traits and compute safety margins
#' @param db.clim formated dataset
#' @param df.traits dataframe of traits associated with eahc species
#' @return dataframe, where for each location, abs/pres of all selected species
#' is mentionned, and associated predictors

get.occurence <- function(db.cont,
                          df.traits,
                          file.path){
  db.cont <- db.cont %>% 
    left_join(df.traits) 
  
  # compute LT50.sp.spring, averaged by species /// IF SPRING FROST NEEDED
  # date.dehardening=0
  # df.fdg.sp <- db.cont %>%
  #   filter(presence==1) %>% 
  #   select(species.binomial,sp.ind,lt50.mean,fdg) %>% 
  #   group_by(species.binomial,sp.ind,lt50.mean) %>% 
  #   summarise(fdg.sp.mean=mean(fdg,na.rm=TRUE),
  #             fdg.sp.sd=sd(fdg,na.rm=TRUE)) %>% 
  #   ungroup() %>% 
  #   mutate(lt50.sp.spring=-5+(30/2)*((5+lt50.mean)/(fdg.sp.mean-date.dehardening)))
  
  
  db.cont <- db.cont %>% 
    # left_join(df.fdg.sp[,c("species.binomial","lt50.sp.spring")],by="species.binomial") %>% 
    mutate(psi_old=psi,
           psi=case_when(psi<(-20000)~(-20000),
                         TRUE~psi),
           psi.100_old=psi.100,
           psi.100=case_when(psi.100<(-20000)~(-20000),
                             TRUE~psi.100)) %>% 
    mutate(hsm=psi-px.mu*1000,
           hsm.100=psi.100-px.mu*1000,
           fsm.winter=frost.winter-lt50.mean,
           # fsm.spring=frost.spring-(lt50.sp.spring),
           species=as.factor(species))
  db.cont$hsm[is.nan(db.cont$hsm)] <- NA
  db.cont$hsm.100[is.nan(db.cont$hsm.100)] <- NA
  db.cont$fsm.winter[is.nan(db.cont$fsm.winter)] <- NA
  # db.cont$fsm.spring[is.nan(db.cont$fsm.spring)] <- NA
  
  
  fwrite(db.cont,file.path)
  
  return(db.cont)
}



#' Define species margins
#' 
#' @description check wether the effect of the occurrence reaches margin limit
#' @param occurence formated dataset
#' @param psi which variable use for psi
#' @param tmin whici variable use for tmin
#' @return dataframe, where for each location, abs/pres of all selected species
#' is mentionned, and associated predictors

get.margins <- function(occurence,
                        psi="psi_cerraday_real",
                        tmin="tmin_cerra"){
  df.quant<-data.frame(species=unique(occurence$species)) |> 
    # left_join(df.species) |> 
    mutate(ymin=NA,
           ymax=NA,
           yabsmin=NA,
           yabsmax=NA,
           quant.psiin=NA,
           quant.psiout=NA,
           quant.tin=NA,
           quant.tout=NA,
           quant.petin=NA,
           quant.petout=NA)
  
  for (sp in df.quant$species){
    
    print(sp)
    path=file.path("data",
                   "chorological_maps_dataset",
                   sp,
                   "shapefiles")
    tryCatch({
      if (file.exists(path)){
        # filter Mauri db with only points of sp
        db.pres <- occurence %>% 
          filter(species==sp) 
        
        #load euforgen distrib
        file.list=grep("plg_clip\\.",
                       list.files(path),
                       value = TRUE)
        file.sp=unique(vapply(strsplit(file.list,"\\."), `[`, 1, FUN.VALUE=character(1)))
        spdistrib=do.call(rbind,lapply(file.sp,function(x)read_sf(dsn=path,x))) |> 
          summarise(geometry = sf::st_union(geometry)) 
        spdistrib=st_crop(spdistrib,
                          xmin=-10,
                          xmax=28.231836,
                          ymin=st_bbox(spdistrib)$ymin[[1]],
                          ymax=st_bbox(spdistrib)$ymax[[1]])
        
        #select all species-points inside euforgen distribution
        db.pres.in <- st_as_sf(db.pres,
                               coords=c("x","y"),
                               crs="epsg:4326") %>% 
          st_join(spdistrib, join = st_within,left=FALSE) %>% # select only points falling in euforgen distrib
          as.data.frame(xy=TRUE) 
        
        # compute lower quantiles of psi and tmin inside distribution
        df.quant[df.quant$species==sp,"quant.psiin"]=quantile(db.pres.in[,psi],
                                                              probs=0.05,
                                                              na.rm=TRUE)[[1]]
        df.quant[df.quant$species==sp,"quant.tin"]=quantile(db.pres.in[,tmin],
                                                            probs=0.05,
                                                            na.rm=TRUE)[[1]]
        df.quant[df.quant$species==sp,"quant.petin"]=quantile(db.pres.in$pet,
                                                              probs=0.95,
                                                              na.rm=TRUE)[[1]]
        
        # select all absences outside euforgen distribution
        db.pres.out <- st_as_sf(db.pres,
                                coords=c("x","y"),
                                crs="epsg:4326") %>% 
          filter(presence==0) |> 
          st_join(spdistrib, join = st_disjoint,left=FALSE) %>% # select only points falling in euforgen distrib
          as.data.frame(xy=TRUE) 
        df.quant[df.quant$species==sp,"quant.psiout"]=quantile(db.pres.out[,psi],
                                                               probs=0.05,
                                                               na.rm=TRUE)[[1]]
        df.quant[df.quant$species==sp,"quant.tout"]=quantile(db.pres.out[,tmin],
                                                             probs=0.05,
                                                             na.rm=TRUE)[[1]]
        df.quant[df.quant$species==sp,"quant.petout"]=quantile(db.pres.out$pet,
                                                               probs=0.95,
                                                               na.rm=TRUE)[[1]]
        
        
        df.quant[df.quant$species==sp,"ymin"]=st_bbox(spdistrib)$ymin[[1]]
        df.quant[df.quant$species==sp,"ymax"]=st_bbox(spdistrib)$ymax[[1]]
        df.quant[df.quant$species==sp,"yabsmin"]=quantile(db.pres[db.pres$presence==0,"y"],
                                                          probs=0.01)[[1]]
        df.quant[df.quant$species==sp,"yabsmax"]=quantile(db.pres[db.pres$presence==0,"y"],
                                                          probs=0.99)[[1]]
        
      } else { # if euforgen distrib do not exist, prevalence default set to 0.1
        print("NA")
      }
    },
    error=function(e){print(paste0("error for ",sp))})
    
  }
  
  df.quant.final<-df.quant |>
    # select(species.binomial,ymin,ymax,yabsmin,yabsmax,
    #        quant.psiin,quant.psiout,quant.tin,quant.tout) |> 
    mutate(hsm.valid.1=yabsmin<ymin,
           fsm.valid.1=yabsmax>ymax,
           hsm.valid.2=quant.psiin>quant.psiout,
           fsm.valid.2=quant.tin>quant.tout,
           hsm.valid.3=quant.petin<quant.petout) 
  
  return(df.quant.final)
}


#' Compute prevalence
#' 
#' @description function that compute for each species in a dataframe its prevalence
#' as the percentage of presence in its presumed distribution
#' @note Distribution from EuForgen
#' @param df.traits df of traits for each species
#' @param db.clim database of pres/abs
#' @return df.loc + extracted variables + computed variables
#' 

get.prevalence <- function(species.list,
                           db.clim){
  df.traits=data.frame(species=species.list)
  # create an empty column for prevalence
  df.traits$prevalence <- NA
  
  # for loop occuring along all species
  
  for (sp in unique(df.traits$species)){
    print(sp)
    # look for euforgen distribution if it exists for the targetted species and
    # compute the prevalelnce
    path=file.path("data",
                   "chorological_maps_dataset",
                   sp,
                   "shapefiles")
    tryCatch(
      {
        if (file.exists(path)){
          db.pres <- db.clim %>%
            filter(species==sp)
          if (length(grep("plg_clip\\.",
                          list.files(path),
                          value = TRUE))!=0){
            #load euforgen distrib
            file.list=grep("plg_clip\\.",
                           list.files(path),
                           value = TRUE)
            file.sp=unique(vapply(strsplit(file.list,"\\."), `[`, 1, FUN.VALUE=character(1)))
          } else if (length(grep("plg\\.",
                                 list.files(path),
                                 value = TRUE))!=0){
            #load euforgen distrib
            file.list=grep("plg\\.",
                           list.files(path),
                           value = TRUE)
            file.sp=unique(vapply(strsplit(file.list,"\\."), `[`, 1, FUN.VALUE=character(1)))
          }
          # filter Mauri db with only points of sp
          spdistrib=do.call(rbind,lapply(file.sp,function(x)read_sf(dsn=path,x))) |> 
            summarise(geometry = sf::st_union(geometry)) 
          spdistrib=st_crop(spdistrib,
                            xmin=-10,
                            xmax=28.231836,
                            ymin=st_bbox(spdistrib)$ymin[[1]],
                            ymax=st_bbox(spdistrib)$ymax[[1]])
          # create spatial points with mauri db 
          db.pres.geo <- st_as_sf(db.pres,
                                  coords=c("x","y"),
                                  crs="epsg:4326") %>% 
            st_join(spdistrib, join = st_within,left=FALSE) %>% # select only points falling in euforgen distrib
            as.data.frame() 
          df.traits[df.traits$species==sp,"prevalence"] <-  sum(db.pres.geo$presence==1)/dim(db.pres.geo)[1]
        } else { # if euforgen distrib do not exist, prevalence default set to 0.1
          df.traits[df.traits$species==sp,"prevalence"] <- 0.1
        }
      },
      error=function(e){print(paste0("error for ",sp))}
    )
  } 
  df.traits$prevalence[is.nan(df.traits$prevalence)] <- mean(df.traits$prevalence,na.rm=TRUE) #because a value is always needed for mod prior
  df.traits$prevalence[is.na(df.traits$prevalence)] <- mean(df.traits$prevalence,na.rm=TRUE)
  return(df.traits[,c("species","prevalence")])
}
