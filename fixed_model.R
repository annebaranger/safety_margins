## libraries
packages.in <- c("stringr","tidyr","dplyr","rstan","foreach","doParallel")
lapply(packages.in,require,character.only=TRUE)

file_occurence="_targets/objects/occurence.binom"
file_df.species="_targets/objects/df.species"
occurence=readRDS(file_occurence)
df.species=readRDS(file_df.species)

species_fit_list=(df.species |> filter(!is.na(px),!is.na(lt50)) |> 
                    filter(species %in% unique(occurence$species)))$species

var.hsm="psi_eraday_real"
var.fsm="tmin_era"
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
  ungroup()

species_correspondance=data.frame(species_num=as.numeric(as.factor(db.clim$species)),
                                  species=db.clim$species) |> 
  unique()

load("mod_random_binom_fixed_alldata.RData")

post<-as.data.frame(t(summary(fit.allsp)$summary))

k_list=post |> select(matches("K_sp")) |> slice(1) |> as.numeric()
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

out=setNames(data.frame(matrix(ncol = 6, nrow = 0)), c("species", "tss", "specificity","sensitivity","thresh_tss","auc"))
for(sp in species_correspondance$species){
  species.num=species_correspondance[species_correspondance$species==sp,
                                     "species_num"]
  K_int=mean(k_list[-species.num])
  db.clim_pred<-db.clim |>
    filter(species==sp) |>
    dplyr::select(species,presence_count,x,y,hsm,fsm,mat,wai) |>
    mutate(presence=as.numeric(presence_count>0),
           pred_sfm=K_int/
             ((1+exp(-post$r_fsm[1]*(fsm-post$t_fsm[1])))*
                (1+exp(-post$r_hsm[1]*(hsm-post$t_hsm[1])))),
           tss=NA,
           thres=NA
    )
  
  observed <- db.clim_pred$presence
  predicted<- db.clim_pred$pred_sfm
  threshold_max=K_int
  thresholds <- seq(0,threshold_max, length.out=100)
  tss_sfm <- sapply(thresholds, calc_tss, observed, predicted)
  out=rbind(out,
            data.frame(species=sp,
                 tss=max(tss_sfm[1,]),
                 specificity=tss_sfm[2,which.max(tss_sfm[1,])],
                 sensitivity=tss_sfm[3,which.max(tss_sfm[1,])],
                 thresh_tss=thresholds[which.max(tss_sfm[1,])],
                 auc=as.numeric(pROC::auc(observed,predicted)))
  )
  
  
}

write.csv(out,file="fit_quality.csv")