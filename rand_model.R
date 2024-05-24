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


# load each species model
file.sandy=list.files("mod.rdata/sandy/")
file.quadri=list.files("mod.rdata/quadri/")[!list.files("mod.rdata/quadri/") %in% file.sandy]
file.list=c(paste0("mod.rdata/sandy/",file.sandy),paste0("mod.rdata/quadri/",file.quadri))
species_issue=c()
file_issue=c()
for(file in file.list){
  load(file)
  if(sum(as.data.frame(summary(fit.allsp)$summary)$Rhat>1.1)>0){
    species_issue=c(species_issue,str_sub(file,1,-7))
    file_issue=c(file_issue,file)
  }
}

file_valid=file.list[!file.list %in% file_issue]
out=setNames(data.frame(matrix(ncol = 8, nrow = 0)), c("species","k_int", "lambda","tss", "specificity","sensitivity","thresh_tss","auc"))

for(file in file_valid){
  load(file)
  post<-as.data.frame(t(summary(fit.allsp)$summary))
  sp<-str_sub(strsplit(file, "/")[[1]][[3]],1,-7)
  db.clim_pred<-db.clim |>
    filter(species==sp) |>
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
  out=rbind(out,
            data.frame(species=sp,
                       k_int=post$K_int[1],
                       lambda=post$lambda[1],
                       tss=max(tss_sfm[1,]),
                       specificity=tss_sfm[2,which.max(tss_sfm[1,])],
                       sensitivity=tss_sfm[3,which.max(tss_sfm[1,])],
                       thresh_tss=thresholds[which.max(tss_sfm[1,])],
                       auc=as.numeric(pROC::auc(observed,predicted)))
  )

  
}

write.csv(out,file="fit_random_quality.csv")
out_random=read.csv("fit_random_quality.csv")
out_fixed=read.csv("fit_quality.csv")
summary(out_random)
summary(out)

t<-out_random |> select(species,tss,specificity,sensitivity,auc)
print(xtable::xtable(t),include.rownames=FALSE)
