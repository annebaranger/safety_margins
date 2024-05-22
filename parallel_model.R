## libraries
packages.in <- c("stringr","tidyr","dplyr","rstan","foreach","doParallel")
lapply(packages.in,require,character.only=TRUE)


## prepare data
file_occurence="occurence.binom"
file_df.species="df.species"
# file_occurence="_targets/objects/occurence.binom"
# file_df.species="_targets/objects/df.species"
occurence=readRDS(file_occurence)
df.species=readRDS(file_df.species)

# species_fit_list=(df.species |> filter(!is.na(px),!is.na(lt50)) |> 
#   filter(species %in% unique(occurence$species)))$species
load("species_to_fit.RData")
species_fit_list=species_to_fit

var.hsm="psi_eraday_real"
var.fsm="tmin_era"
folder.out="mod.rdata"
if(!dir.exists(folder.out)){dir.create(folder.out)}
  # filter(species %in% species.select) |>
# species.select=sample(unique(db.clim$species),20)
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
  group_by(species) |>
  sample_frac(0.6) |> #!!
  ungroup()


## functions
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


fit_random_sp<-function(db.clim,
                        sp.excl
                        ){
  db.sp <- db.clim |> filter(species!=sp.excl)
  data.list<-list(N=dim(db.sp)[1],
                  S=nlevels(as.factor(db.sp$species)),
                  max_draws=max(db.sp$n),
                  presence=db.sp$presence_count,
                  draw=db.sp$n,
                  species=as.numeric(as.factor(db.sp$species)),
                  fsm=db.sp$fsm,
                  hsm=db.sp$hsm)
  print(paste0("Launching model for species : ",sp.excl))

  fit.allsp <- stan(file = "glm_log_all_betapareto.stan",
                    data=data.list,
                    iter=1000, #!!
                    # core=3,#!!
                    chains=3,
                    include=FALSE,
                    pars=c("proba","K_vect"))

  file_path=file.path(folder.out,paste0(sp.excl,".rdata"))
  save(fit.allsp,file=file_path)
  # load("mod_random_binom_betapareto_alldata.RData")
  
  print(paste0("Post-process for species : ",sp.excl))
  

  post<-as.data.frame(t(summary(fit.allsp)$summary)) |>
    select(!matches("K_sp"))


  db.clim_pred<-db.clim |>
    filter(species==sp.excl) |>
    dplyr::select(species,presence_count,x,y,hsm,fsm,mat,wai) |>
    mutate(presence=as.numeric(presence_count>0),
           pred_sfm=post$K_int[1]/
             ((1+exp(-post$r_fsm[1]*(fsm-post$t_fsm[1])))*
                (1+exp(-post$r_hsm[1]*(hsm-post$t_hsm[1])))),
           tss=NA,
           thres=NA
    )

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
                 rhat=post$lp__[10])
  return(out)
}
# fit_random_sp(db.clim,sp.excl = "Abies alba")

#Setup backend to use many processors
print("launching parallelisation")
totalCores = detectCores()

#Leave one core to avoid overload your computer
cluster <- makeCluster(min(totalCores[1]-2,38)) 
registerDoParallel(cluster)
clusterExport(cluster, c("db.clim", "fit_random_sp", "df.species","calc_tss","folder.out"))

out <- foreach(i = species_fit_list, .combine=rbind,.packages = packages.in) %dopar% {
  library(rstan)
  library(Rcpp)
  fit_random_sp(db.clim=db.clim,sp.excl=i)
  # print(i)
}

write.csv(out,file=file.path(folder.out,"out.csv"))
parallel::stopCluster(cluster)

