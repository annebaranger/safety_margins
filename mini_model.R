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

species_fit_list=(df.species |> filter(!is.na(px),!is.na(lt50)) |> 
                    filter(species %in% unique(occurence$species)))$species

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
sp.excl="Abies alba"
db.sp <- db.clim |> filter(species!=sp.excl)

data.list<-list(N=dim(db.sp)[1],
                S=nlevels(as.factor(db.sp$species)),
                max_draws=max(db.sp$n),
                presence=db.sp$presence_count,
                draw=db.sp$n,
                species=as.numeric(as.factor(db.sp$species)),
                fsm=db.sp$fsm,
                hsm=db.sp$hsm)
