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
                       dws=FALSE,
                       file.path){
  ## Load safe ty margin indicators
  clim.list=lapply(clim.list,rast)
  df.clim=lapply(seq_along(clim.list),function(x) {
    df=terra::extract(clim.list[[x]],data.frame(x=db.mauri$x,y=db.mauri$y))[-1]
    colnames(df)=names(clim.list)[x]
    return(df)}
  )
  
  db.mauri<-cbind(db.mauri,
                  as.data.frame(df.clim))
  
  db.mauri$presence[is.na(db.mauri$presence)] <- 0
  
  
  if(dws){
    var_names=c("cell",names(clim.list)[names(clim.list)!="cell"],"species")
    db.mauri<-db.mauri |> 
      group_by_at(var_names) |> 
      summarize(presence_count=sum(presence,na.rm = TRUE),
                n=n(),
                x=mean(x),
                y=mean(y)) |> 
      ungroup()
  }
  
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

#' Compute niche caracteristics
#' 
#' @description compute mean lat/long and continentality of species niche
#' @note Continentality from https://doi.org/10.1038/s41597-020-0464-0
#' @param df.traits df of traits for each species
#' @param db.clim database of pres/abs
#' @return df.loc + extracted variables + computed variables
#' 

get.niche <- function(species.list,
                      db.clim,
                      psi="psi_cerraday_real",
                      tmin="tmin_cerra"){
  df.traits=data.frame(species=species.list)
  
  df.niche <- db.clim %>% 
    filter(presence==1) %>% 
    group_by(species) %>% 
    summarise(lat.mean=mean(y),
              lat.sd=sd(y),
              long.mean=mean(x),
              long.sd=sd(x),
              lat.q05=quantile(y,prob=0.05)[[1]],
              lat.q95=quantile(y,prob=0.95)[[1]],
              psi.q05=quantile(eval(parse(text=psi)),prob=0.05,na.rm=TRUE)[[1]],
              psi.q95=quantile(eval(parse(text=psi)),prob=0.95,na.rm=TRUE)[[1]],
              tmin.q05=quantile(eval(parse(text=tmin)),prob=0.05,na.rm=TRUE)[[1]],
              tmin.q95=quantile(eval(parse(text=tmin)),prob=0.95,na.rm=TRUE)[[1]]
    )
  
  
  # rast.cont=as.data.frame(mean(rast("data/jci_year.nc")),xy=TRUE) %>% 
  #   mutate(z=x,
  #          x=y,
  #          y=z) %>% 
  #   select(-z) %>% 
  #   rast(crs="epsg:4326")
  # df.traits$jci=NA
  # for (sp in df.traits$species.binomial){
  #   print(sp)
  #   db.pres <- db.clim %>% 
  #     filter(species.binomial==sp) %>% 
  #     filter(presence==1)
  #   db.pres <- cbind(db.pres,
  #                    jci=extract(rast.cont,db.pres[,c("x","y")])[["mean"]])
  #   jci=mean(db.pres$jci,na.rm=TRUE)
  #   df.traits[df.traits$species.binomial==sp,"jci"]=jci
  # }
  
  # df.overlap <- db.clim %>%
  #   filter(presence==1) %>%
  #   group_by(node) %>%
  #   mutate(nb.sp=sum(presence==1),
  #          overlap.hsm=sum(hsm>0),
  #          overlap.fsm=sum(fsm.winter>0)) %>%
  #   ungroup() %>%
  #   group_by(species.binomial) %>%
  #   summarise(overlap=mean(nb.sp),
  #             overlap.hsm=mean(overlap.hsm,na.rm=TRUE),
  #             overlap.fsm=mean(overlap.fsm,na.rm=TRUE))
  
  df.traits <- df.traits %>%
    left_join(df.niche,by="species") #%>%
  # left_join(df.overlap,by="species.binomial")
  return(df.traits[,c("species","lat.mean","lat.sd","long.mean","long.sd",
                      "lat.q05","lat.q95","psi.q05","psi.q95","tmin.q05","tmin.q95")])#,"jci","overlap","overlap.hsm","overlap.fsm"
}


#' Get shade tol
#' 
#' @description compute mean lat/long and continentality of species niche
#' @note Continentality from https://doi.org/10.1038/s41597-020-0464-0
#' @param df.traits df of traits for each species
#' @return df.loc + extracted variables + computed variables
#' 

get.shadetol <- function(species.list,
                         db.clim){
  df.traits=data.frame(species=species.list)
  df.shadetol=read.csv2("data/Species traits/data_Niinemets&Valladares_2006.csv")
  df.traits=df.traits %>% 
    left_join(df.shadetol,by=c("species"="Species")) |> 
    mutate(across(c(shade_tolerance.mean,drought_tolerance.mean,waterlogging_tolerance.mean),
                  as.numeric))
  
  df.overlapshade=db.clim |>
    filter(presence==1) |> 
    select(node,species,presence) |> 
    left_join(df.traits[,c("species","shade_tolerance.mean","drought_tolerance.mean")], by="species") |> 
    filter(!is.na(shade_tolerance.mean)&
             !is.na(drought_tolerance.mean)) |> 
    group_by(node) |> 
    mutate(overlap_plot_shade=colSums(outer(shade_tolerance.mean,shade_tolerance.mean,">")),
           overlap_plot_drought=colSums(outer(drought_tolerance.mean,drought_tolerance.mean,">"))) |> #shade_overlap(shade_tolerance.mean,presence) 
    ungroup() |> 
    group_by(species) |> 
    summarise(overlap_shade=mean(overlap_plot_shade),
              overlap_drought=mean(overlap_plot_drought))
  
  df.traits=df.traits %>% 
    left_join(df.overlapshade,by=c("species"))
  return(df.traits[,c("species","shade_tolerance.mean","drought_tolerance.mean","waterlogging_tolerance.mean","overlap_shade","overlap_drought")]) 
}

#' Get traits and caract
#' 
#' @description compute niche traits and caract
#' @param df.traits df of traits for each species
#' @return df.loc + extracted variables + computed variables
#' 

get.species <- function(species.list,
                        df.preval,
                        df.shadetol,
                        df.niche,
                        df.traits,
                        file.output){
  df.species <- data.frame(species=species.list) %>% 
    left_join(df.shadetol,by="species")%>% 
    left_join(df.niche,by="species")%>% 
    left_join(df.traits,by="species") |> 
    left_join(df.preval,by="species")
  fwrite(df.species,file=file.output)
  return(df.species) 
}


#' Fit random model with cross validation
#' 
#'@description fit random model for each species
#'@param db.clim
#'@param df.traits
#'@param sp.excl
#'@return file of model fit

fit_random_sp<-function(occurence,
                        df.species,
                        var.hsm="psi_eraday_real",
                        var.fsm="tmin_era",
                        sp.excl,
                        folder.out="mod.rdata"){
  print(sp.excl)
  if(!dir.exists(folder.out)){dir.create(folder.out)}
  
 
  
  db.clim=occurence %>% 
    ungroup() |> 
    left_join(df.species) |> 
    mutate(hsm:=(!!sym(var.hsm)/1000)-px,
           fsm:=!!sym(var.fsm)-lt50)  |> 
    filter(!is.na(hsm)) |> 
    filter(!is.na(fsm)) |> 
    filter(hsm>quantile(hsm,prob=0.01)) |> 
    filter(fsm>quantile(hsm,prob=0.01)) |> 
    filter(!is.na(wai)) |> 
    filter(!is.na(mat)) |> 
    filter(species!=sp.excl) 
  species.select=sample(unique(db.clim$species),20)
  db.clim=db.clim |> 
    filter(species %in% species.select) |> 
    group_by(species) |> 
    sample_n(20000) |> 
    ungroup()
  
  data.list<-list(N=dim(db.clim)[1],
                  S=nlevels(as.factor(db.clim$species)),
                  max=max(db.clim$n),
                  presence=db.clim$presence_count,
                  draw=db.clim$n,
                  species=as.numeric(as.factor(db.clim$species)),
                  fsm=db.clim$fsm,
                  hsm=db.clim$hsm)
 
  fit.allsp <- stan(file = "glm_log_all.stan",
                    data=data.list,
                    # init=init,
                    iter=500,
                    chains=3,
                    core=3,
                    include=FALSE,
                    pars=c("proba","K_vect"))
  
  file_path=file.path(folder.out,paste0(sp.excl,".rdata"))
  save(file_path)
  
  post<-as.data.frame(t(summary(fit.allsp)$summary)) |> 
    select(!matches("K_sp"))
  
  
  db.clim_pred<-occurence %>% 
    ungroup() |> 
    left_join(df.species) |> 
    mutate(hsm:=(!!sym(var.hsm)/1000)-px,
           fsm:=!!sym(var.fsm)-lt50)  |> 
    filter(!is.na(hsm)) |> 
    filter(!is.na(fsm)) |> 
    filter(hsm>(-10000)) |> 
    filter(!is.na(wai)) |> 
    filter(!is.na(mat)) |> 
    filter(species==sp.excl) |> 
    # filter(species.binomial==sp) |> 
    select(species,presence_count,x,y,hsm,fsm,mat,wai) |> 
    mutate(presence=as.numeric(presence_count>0),
           pred_sfm=post$K_int[1]/
             ((1+exp(-post$r_fsm[1]*(fsm-post$t_fsm[1])))*
                (1+exp(-post$r_hsm[1]*(hsm-post$t_hsm[1])))),
           tss=NA,
           thres=NA
    ) 
  
  calc_tss <- function(threshold, observed, predicted_probs) {
    observed=as.factor(observed)
    predicted <- ifelse(predicted_probs > threshold, 1, 0)
    predicted=factor(predicted,levels=c(0,1))
    conf_matrix <- table(observed, predicted)
    TP <- conf_matrix[2, 2]
    FN <- conf_matrix[2, 1]
    TN <- conf_matrix[1, 1]
    FP <- conf_matrix[1, 2]
    sensitivity <- TP / (TP + FN)
    specificity <- TN / (TN + FP)
    TSS <- sensitivity + specificity - 1
    return(c(TSS,specificity,sensitivity))
  }
  
  observed <- db.clim_pred$presence
  predicted<- db.clim_pred$pred_sfm
  threshold_max=post$K_int[1]
  thresholds <- seq(0,threshold_max, length.out=100)
  tss_sfm <- sapply(thresholds, calc_tss, observed, predicted)
  out=data.frame(species=sp.excl,
                 file=file_path,
                 tss=max(tss_sfm[1,]),
                 specificity=tss_sfm[2,which.max(tss_sfm[1,])],
                 sensitivity=tss_sfm[3,which.max(tss_sfm[1,])],
                 thresh_tss=thresholds[which.max(tss_sfm[1,])],
                 rhat=post$rhat[1])
  
  return(out)
  
}
