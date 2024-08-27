## packages
packages.in <- c("targets","stringr","tidyr","dplyr","rstan","ggplot2","sf",
                 "latex2exp","ggpubr","ggsignif","terra")
lapply(packages.in,require,character.only=TRUE)

## loads
tar_load(df.species) 
db.clim<-tar_read(occurence) |> 
  left_join(df.species) |> 
  mutate(hsm=(psi_eraday_real/1000)-px,
         fsm=tmin_era-lt50) 

tar_load(europe)
# df.species <- df.species |>
#   left_join(
#     db.clim |> filter(presence==1) |> 
#       group_by(species) |> 
#       summarise(n=n()))

# tar_load(traits)

tar_load(species_fit_list)
tar_load(safmarg_sp_select)
# tar_load(mod.select.safmarg_beta)
tar_load(margin.limit)


#%%%%%%%%%%%%%%
#### fig 1 ####
#%%%%%%%%%%%%%%
predictor=c("tmin_era","psi_eraday_real","wai","pet","map","mat","prpet","sgdd")

# compute correlations
db.pred<-db.clim |> 
  select(x,y,map,pet,psi_eraday_real,mat,tmin_era) |> 
  distinct() |> 
  drop_na()
cor.test(db.pred$pet,db.pred$psi_eraday_real)
cor.test(db.pred$map,db.pred$psi_eraday_real)
cor.test(db.pred$mat,db.pred$tmin_era)

# get climatic data and traits by species 
db.clim.lm <- db.clim |>  
  filter(presence==1) |> 
  filter(psi_eraday_real>(-15000)) |>
  left_join(margin.limit) |> 
  mutate(prpet=map-pet) |> 
  pivot_longer(cols=all_of(predictor)) |> 
  group_by(name) |> 
  mutate(value=scale(value)) |> 
  ungroup() |> 
  group_by(species,hsm.valid.3,fsm.valid.2,name,lt50,px) |> 
  summarise(quant05=quantile(value,probs=0.05,na.rm=TRUE)[[1]],
            quant5=quantile(value,probs=0.5,na.rm=TRUE)[[1]],
            quant95=quantile(value,probs=0.95,na.rm=TRUE)[[1]]) |> 
  ungroup() |> 
  mutate(lt50=scale(lt50)[,1],
         px=scale(px)[,1]) |> 
  pivot_longer(cols=c("quant05","quant5","quant95"),
               names_to="quant_name",
               values_to = "quant_val") 

# raw dataframe
df.lm.traits<-data.frame(pred=rep(predictor,2)) |> 
  crossing(trait=c("lt50","px"),
           quant=c("quant05","quant95"),
           mean=NA,
           sd=NA,
           pval=NA,
           r2=NA,
           conf.up=NA,
           conf.low=NA,
           n=NA)


for(pred in predictor){
  print(pred)
  for(quant in c("quant05","quant95")){
    print(quant)
    
    # test for p50
    db.clim.pred<-db.clim.lm |> 
      filter(name==pred) |> 
      filter(quant_name==quant)# |> 
    # filter(hsm.valid.3==TRUE) 
    # filter(mod%in%c("2sm","hsm"))
    summary=as.data.frame(summary(lm(quant_val~px,db.clim.pred))$coefficients)
    conf=confint(lm(formula = quant_val ~ px, data = db.clim.pred))
    print(summary(lm(quant_val~px,db.clim.pred)))
    df.lm.traits[df.lm.traits$pred==pred&
                   df.lm.traits$quant==quant&
                   df.lm.traits$trait=="px",
                 c("mean","sd","pval","r2","conf.up","conf.low","n")]=
      c(summary["px",c("Estimate","Std. Error","Pr(>|t|)")],
        summary(lm(quant_val~px,db.clim.pred))$r.squared,
        conf["px",1],
        conf["px",2],
        dim(db.clim.pred)[1])
    
    
    #test for lt50
    db.clim.pred<-db.clim.lm |> 
      filter(name==pred) |> 
      filter(quant_name==quant)# |> 
    # filter(fsm.valid.2==TRUE)
    # filter(mod%in%c("2sm","fsm"))
    summary=as.data.frame(summary(lm(quant_val~lt50,db.clim.pred))$coefficients)
    conf=confint(lm(formula = quant_val ~ lt50, data = db.clim.pred))
    print(summary(lm(quant_val~lt50,db.clim.pred)))
    df.lm.traits[df.lm.traits$pred==pred&
                   df.lm.traits$quant==quant&
                   df.lm.traits$trait=="lt50",
                 c("mean","sd","pval","r2","conf.up","conf.low","n")]=
      c(summary["lt50",c("Estimate","Std. Error","Pr(>|t|)")],
        summary(lm(quant_val~lt50,db.clim.pred))$r.squared,
        conf["lt50",1],
        conf["lt50",2],
        dim(db.clim.pred)[1])
    rm(db.clim.pred)
  }
}


df.lm.traits<-df.lm.traits |> 
  mutate(signif1=pval<0.001,
         signif2=pval>0.001&pval<0.01,
         signif3=pval>0.01&pval<0.05,
         text1=ifelse(signif3,"*",
                      ifelse(signif2,"**",
                             ifelse(signif1,"***",""))),
         text2=ifelse(signif3, paste0("R2=",round(r2, digits = 2)),
                      ifelse(signif2,paste0("R2=",round(r2,digits=2)),
                             ifelse(signif1,paste0("R2=",round(r2,digits=2)),""))))

df.lm.traits |> 
  # filter the relevant quantile for each variable (ie stressful edge)
  filter((pred%in% c("sgdd","tmin_era","psi_eraday_real","wai","prpet","map","mat") & quant=="quant05")|
           (pred=="pet"& quant=="quant95")) |>
  # filter the relevant trait for each variable
  filter((trait=="lt50"&(pred %in% c("sgdd","tmin_era","mat")))|
           (trait=="px"&(pred %in% c("wai","psi_eraday_real","prpet","map","pet")))) |> 
  # rename traits
  mutate(trait=case_when(trait=="px"~"$\\Psi_{crit}$",
                         trait=="lt50"~"$LT_{50}$")) |> 
  # filter variables
  filter(pred %in% c("mat","tmin_era","psi_eraday_real","pet","map")) |> 
  # rename
  mutate(pred=case_when(pred=="tmin_era"~"$T_{min}$",
                        pred=="psi_eraday_real"~"$\\Psi_{min}$",
                        TRUE~pred)) |> 
  mutate(quant=forcats::fct_recode(quant, "5%"="quant05", "95%"="quant95")) |> 
  ggplot(aes(x=mean,y=pred,color=quant))+
  geom_point(size=1.5,alpha=0.6,position=position_dodge(width=0.5))+ #,position=position_dodge(width=0.5)
  geom_pointrange(aes(xmin=conf.low,xmax=conf.up),position=position_dodge(width=0.5))+ #,position=position_dodge(width=0.5)
  geom_text(aes(label = text1,
                x=mean,
                y=pred
  ),
  # position = position_dodge(9),
  size=15/.pt,
  show.legend = FALSE,
  vjust=-0.25)+
  geom_text(aes(label = text2,
                x=mean,
                y=pred
  ),
  # position = position_dodge(9),
  size=8/.pt,
  show.legend = FALSE,
  vjust=2)+
  geom_vline(xintercept = 0,linetype="dotted")+
  scale_color_manual(values=c("darkred","lightblue"))+ 
  scale_y_discrete(label=TeX)+
  facet_wrap(~trait,
             labeller=as_labeller(TeX,
                                  default = label_parsed),
             scales = "free_y")+
  theme_minimal()+
  theme(axis.title= element_blank(),
        axis.ticks.x = element_line(linewidth = 0.6),
        axis.line.x = element_line(linewidth = 0.6),
        panel.grid = element_blank(),
        legend.position = "top",
        axis.text.y = element_text(hjust=0.5),
        text=element_text(size=11))+
  xlim(-1,1.1)+
  labs(color="Quantiles")->fig1


#%%%%%%%%%%%%%%
#### fig 2 ####
#%%%%%%%%%%%%%%
hsm.95=quantile(db.clim$hsm,prob=0.95,na.rm=TRUE)[[1]]
hsm.05=quantile(db.clim$hsm,prob=0.05,na.rm=TRUE)[[1]]
fsm.95=quantile(db.clim$fsm,prob=0.95,na.rm=TRUE)[[1]]
fsm.05=quantile(db.clim$fsm,prob=0.05,na.rm=TRUE)[[1]]

df.range<-db.clim |> 
  filter(presence==1) |> 
  filter(psi_eraday_real>(-10000)) |> 
  group_by(species) |> 
  summarise(psi_range=max(psi_eraday_real)/1000-min(psi_eraday_real)/1000,
            t_range=max(tmin_era,na.rm=TRUE)-min(tmin_era,na.rm=TRUE))

#first panel
c<-safmarg_sp_select %>%
  left_join(df.species) |>
  crossing(data.frame(xsm_val=c(seq(fsm.05,
                                    fsm.95,
                                    length.out=100),
                                seq(hsm.05,
                                    hsm.95,
                                    length.out=100)),
                      xsm_name=c(rep("Frost safety margins (°C)",100),
                                 rep("Hydraulic safety margins (MPa)",100)))) |> 
  filter((mod %in% c("2var","hsm") & xsm_name=="Hydraulic safety margins (MPa)")|
           mod %in% c("2var","fsm") & xsm_name=="Frost safety margins (°C)") |> 
  mutate(pred=case_when(mod=="2var"&xsm_name=="Hydraulic safety margins (MPa)"~k_int/((1+exp(-r_fsm*(fsm.95-t_fsm)))*
                                                                                        (1+exp(-r_hsm*(xsm_val-t_hsm)))),
                        mod=="2var"&xsm_name=="Frost safety margins (°C)"~k_int/((1+exp(-r_hsm*(hsm.95-t_hsm)))*
                                                                                   (1+exp(-r_fsm*(xsm_val-t_fsm)))),
                        mod=="hsm"~k_int/(1+exp(-r_hsm*(xsm_val-t_hsm))),
                        mod=="fsm"~k_int/(1+exp(-r_fsm*(xsm_val-t_fsm)))
  )
  ) |>
  ggplot(aes(xsm_val,pred,color=taxa,linesize=species,linetype=mod))+ #color=prevalence,
  geom_vline(xintercept = 0,linetype="dashed")+
  geom_line(alpha=0.6,size=0.6)+
  scale_color_manual(values=c("chocolate2","forestgreen"))+
  # scale_color_gradientn(colours = viridis(15,direction=-1),
  #                       breaks=c(0,0.1,0.2,0.3,0.4,0.5,0.6))+
  guides(alpha="none")+
  labs(color="Species prevalence")+
  theme_minimal()+
  theme(axis.title= element_blank(),
        axis.ticks.x = element_line(linewidth = 0.6),
        axis.line = element_line(linewidth = 0.6),
        legend.position = "none",
        # panel.grid = element_blank(),
        # panel.grid.major = element_line(color = "lightgrey",
        #                                   size =0.02),
        axis.text.y = element_text(hjust=0.5),
        strip.placement = "outside",
        strip.text = element_text(size=11,vjust=-1.1),
        text=element_text(size=11))+
  facet_wrap(~xsm_name,scales="free",
             strip.position = "bottom")+
  NULL

#second panel
df.n<- safmarg_sp_select %>% 
  left_join(df.species) |> 
  rename(HSM=inflex_hsm,
         FSM=inflex_fsm) %>% 
  pivot_longer(col=c("HSM","FSM")) %>% 
  filter(mod!="none") %>% 
  filter(!(mod=="fsm"&name=="HSM")&
           !(mod=="hsm"&name=="FSM")) %>% 
  group_by(taxa,name) |> 
  summarise(n=n(),
            value=quantile(value,probs=0.5)[[1]]) |> 
  ungroup() |> 
  mutate(value=case_when(name=="FSM"~95,
                         TRUE~value))
safmarg_sp_select %>%
  left_join(df.species) |> 
  mutate(HSM=inflex_hsm,
         FSM=inflex_fsm) %>% 
  pivot_longer(col=c("HSM","FSM")) %>% 
  filter(mod!="none") %>% 
  filter(!(mod=="fsm"&name=="HSM")&
           !(mod=="hsm"&name=="FSM")) %>% 
  ggplot(aes(name,value,fill=taxa))+
  geom_boxplot()+
  geom_text(data=df.n,
            mapping=aes(x=name,y=value,label=n),
            position = position_dodge(width=0.75),
            vjust=2.5,
            size=9/.pt)+
  geom_signif(comparisons = list(c("HSM","FSM")),
              map_signif_level = TRUE)+
  # geom_signif(comparisons = list(c("angiosperm","gymnosperm")),
  #             map_signif_level = TRUE)+
  # stat_compare_means(aes(label = ..p.signif..),
  #                    method = "t.test",
  #                    label.x = 1.5, label.y = 60)+
  scale_fill_manual(values=c("chocolate2","forestgreen"))+
  theme_minimal()+
  theme(axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_line(linewidth = 0.6),
        legend.title = element_blank(),
        # panel.grid = element_blank(),
        # panel.grid.major = element_line(color = "lightgrey",
        #                                   size =0.02),
        axis.text.y = element_text(hjust=0.5),
        text=element_text(size=11))+
  labs(x="Inflexion index")-> b
#%%%%%%%%%%%%%%%%%%%
#third panel
df.n<- safmarg_sp_select %>% 
  left_join(df.species) |>
  left_join(df.range) |> 
  dplyr::select(species,px,lt50,t_hsm,t_fsm,mod,psi_range,t_range,taxa) %>% 
  mutate(t_range=case_when(is.na(t_range)~mean(df.range$t_range,na.rm=TRUE),
                           TRUE~t_range),
         HSM=t_hsm/psi_range,
         FSM=t_fsm/t_range) |> 
  dplyr::select(taxa,species,mod,HSM,FSM) %>% 
  pivot_longer(cols=c("HSM","FSM"),names_to = "threshold",values_to = "traits_val") %>% 
  filter((mod%in%c("hsm","2var")&threshold=="HSM")|
           (mod%in%c("fsm","2var")&threshold=="FSM")) %>% 
  group_by(taxa,threshold) |> 
  summarise(n=n(),
            value=quantile(traits_val,probs=0.5)[[1]]) |> 
  ungroup()
safmarg_sp_select %>% 
  left_join(df.species) |>
  left_join(df.range) |> 
  dplyr::select(species,px,lt50,t_hsm,t_fsm,mod,psi_range,t_range,taxa) %>% 
  mutate(t_range=case_when(is.na(t_range)~mean(df.range$t_range,na.rm=TRUE),
                           TRUE~t_range),
         `HSM`=t_hsm/psi_range,
         FSM=t_fsm/t_range) |> 
  pivot_longer(cols=c("px","lt50"),names_to = "traits",values_to = "traits_val") %>% 
  pivot_longer(cols=c("HSM","FSM"),names_to="threshold",values_to="threshold_val") %>% 
  filter((mod%in%c("hsm","2var")&traits=="px"&threshold=="HSM")|
           (mod%in%c("fsm","2var")&traits=="lt50"&threshold=="FSM")) %>% 
  group_by(threshold) |> 
  mutate(n=n()) |> 
  ungroup() |> 
  ggplot(aes(threshold,threshold_val,fill=taxa))+
  geom_boxplot()+
  geom_text(data=df.n,
            mapping=aes(x=threshold,y=value,label=n),
            position = position_dodge(width=0.75),
            vjust=-0.5,
            size=8/.pt)+
  geom_signif(comparisons = list(c("HSM","FSM")),
              map_signif_level = TRUE)+
  # geom_text(data=df.n,mapping=aes(x=threshold,y=threshold_val,label=n))+
  # stat_compare_means(aes(label = ..p.signif..),
  #                    method = "t.test",
  #                    label.x = 1.5, label.y = 0.35)+
  scale_fill_manual(values=c("chocolate2","forestgreen"))+
  theme_minimal()+
  theme(axis.title.y= element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_line(linewidth = 0.6),
        # panel.grid = element_blank(),
        # panel.grid.major = element_line(color = "lightgrey",
        #                                   size =0.02),
        axis.text.y = element_text(hjust=0.5),
        strip.placement = "outside",
        text=element_text(size=11))+
  labs(x="Rescaled thresholds")+
  NULL-> a

cowplot::plot_grid(
  c,
  cowplot::plot_grid(get_legend(b),
                     cowplot::plot_grid(a+theme(legend.position = "none"),
                                        b+theme(legend.position = "none"),
                                        nrow=1),
                     nrow=2,
                     rel_heights = c(0.2,1),
                     rel_widths = c(0.2,1),
                     axis="l"
  ),
  labels=c("A","B"),
  label_size = 12,
  hjust=c(-0.2,0.51),
  rel_widths = c(1.1,0.6)
)->fig2


#%%%%%%%%%%%%%%%%
#### Sup mat ####
#%%%%%%%%%%%%%%%%

### S1
tmin=rast(tar_read(tmin2m_era))
names(tmin)="tmin"
as.data.frame(tmin,
              xy=TRUE) |> 
  mutate(tmin=cut(tmin,
                  breaks=c(-Inf, -50,-40, -30,-20,-15,-10,-5,0,Inf),
                  labels=c("<-50", "-50<Tmin<-40","-40/-30","-30/-20","-20/-15","-15/-10","-10/-5",
                           "-5/0",">0"))) %>%
  ggplot() +
  geom_tile(aes(x=x,y=y,fill=tmin))+
  geom_sf(data=europe,
          fill=NA)+
  coord_sf()+
  theme_minimal() +
  labs(fill="Tmin")+
  theme(axis.title=element_blank(),
        axis.text = element_blank(),
        legend.key.size = unit(0.5,"cm"),
        legend.text = element_text(size=8),
        legend.title = element_text(size=9))+
  scale_fill_brewer(palette="RdYlBu")->figS1

### S2

theta=seq(0,1,by=0.00001)
hydrau_pars=data.frame(hor=as.numeric(),alpha=as.numeric(),m=as.numeric(),n=as.numeric(),thr=as.numeric(),ths=as.numeric())
for(i in 1:7){
  print(i)
  alpha=global(rast(paste0("data/EU_SoilHydroGrids_1km/MRC_alp_sl",i,".tif")),'mean',na.rm=TRUE)[[1]]
  m=global(rast(paste0("data/EU_SoilHydroGrids_1km/MRC_m_sl",i,".tif")),'mean',na.rm=TRUE)[[1]]
  n=global(rast(paste0("data/EU_SoilHydroGrids_1km/MRC_n_sl",i,".tif")),'mean',na.rm=TRUE)[[1]]
  thr=global(rast(paste0("data/EU_SoilHydroGrids_1km/MRC_thr_sl",i,".tif")),'mean',na.rm=TRUE)[[1]]
  ths=global(rast(paste0("data/EU_SoilHydroGrids_1km/MRC_ths_sl",i,".tif")),'mean',na.rm=TRUE)[[1]]
  hydrau_pars=rbind(hydrau_pars,list(hor=i,alpha=alpha,m=m,n=n,thr=thr,ths=ths))
}
merge(hydrau_pars,
      as.data.frame(theta)) |> 
  mutate(REW=(theta-thr*10^(-4))/((ths-thr)*10^(-4)),
         psi=case_when(theta<ths~(((1/REW)^(1/(m*10^(-4)))-1)^(1/(n*10^(-4)))*(1/(alpha*10^(-4)))*9.78*10^(-2))/1000,
                       theta>ths~ths*10^(-4)),
         Horizon=as.factor(hor)) |> 
  filter(psi>0.00002) |> 
  ggplot(aes(theta,psi,color=Horizon))+
  geom_line()+
  scale_y_log10()+
  scale_color_brewer(palette="Spectral",direction=-1)+
  labs(x=TeX("$\\theta$ (%)"),
       y=TeX("$\\Psi$ (MPa)"))+
  theme_minimal()->figS2

### S3
psi=aggregate(c(rast(tar_read(psi_eradaycdo_real_beta)[,c("x","y","psi")]),
                rast(tar_read(psi_eradaycdo_real)[,c("x","y","psi")])),
              10,
              method="mean",
              na.rm=TRUE)
names(psi)=c("psi_beta","psi_real")


as.data.frame(psi) |> 
  mutate(ssy=((psi_beta-psi_real)/abs(psi_beta))) |> 
  summarise(ssy=median(ssy,na.rm=TRUE),
            cor=cor(psi_beta,psi_real, use="pairwise.complete.obs"),
            rmse=exp(sqrt(mean((log(-psi_beta)-log(-psi_real))^2,na.rm=TRUE))))

t_test=as.data.frame(psi) |> 
  mutate(ssy=((psi_beta-psi_real)/abs(psi_beta))) |> 
  drop_na() |> 
  summarise(t.test=t.test(psi_beta,psi_real)$p.value)


as.data.frame(psi,xy=TRUE)|>
  mutate(dif=cut((psi_real-psi_beta)/1000,
                 breaks=c(-Inf, -2, -1, -0.5, -0.25, 0, 0.25,0.5,1,2,Inf))#,
         # psi_beta=cut(psi_beta,
         #       breaks=c(-Inf, -15000, -10000, -5000,-3000, -2000, -1000, -500, 0,Inf),
         #               labels=c("<-15MPa","-15<HSM<-10MPa","-10<HSM<-5MPa",
         #                        "-5<psi<-3MPa","-3/-2","-2/-1","-1/-0.5","-0.5/0", ">0MPa")),
  ) |> 
  # pivot_longer(cols=c("psi_real","psi_beta")) |> 
  ggplot()+
  geom_tile(aes(x=x,y=y,fill=dif))+
  geom_sf(data=europe,
          fill=NA)+
  coord_sf()+
  theme_minimal() +
  theme(axis.title=element_blank(),
        axis.text = element_blank(),
        legend.key.size = unit(0.5,"cm"),
        legend.text = element_text(size=10),
        legend.title = element_text(size=11))+
  scale_fill_brewer(palette="RdYlBu")+
  labs(fill=unname(TeX("$\\Psi_{min}-\\Psi_{new}(MPa)$")))-> figS3


### S4
psi=aggregate(c(rast(tar_read(sensitivity_2_0.966_output.psi_era_day_real_2_0.966)[,c("x","y","psi")]),
                rast(tar_read(sensitivity_3_0.966_output.psi_era_day_real_3_0.966)[,c("x","y","psi")]),
                rast(tar_read(sensitivity_4_0.966_output.psi_era_day_real_4_0.966)[,c("x","y","psi")]),
                rast(tar_read(sensitivity_5_0.966_output.psi_era_day_real_5_0.966)[,c("x","y","psi")]),
                rast(tar_read(sensitivity_6_0.966_output.psi_era_day_real_6_0.966)[,c("x","y","psi")]),
                rast(tar_read(sensitivity_7_0.966_output.psi_era_day_real_7_0.966)[,c("x","y","psi")]),
                rast(tar_read(sensitivity_8_0.966_output.psi_era_day_real_8_0.966)[,c("x","y","psi")])
                ),
              10,
              method="mean",
              na.rm=TRUE)
names(psi)=c("psi_2","psi_3","psi_4","psi_5","psi_6","psi_7","psi_8")

## summary stats
as.data.frame(psi,xy=TRUE)|>
  pivot_longer(cols=c("psi_2","psi_3","psi_4","psi_6","psi_7","psi_8")) |> 
  mutate(ssy=((psi_5-value)/abs(psi_5))) |> 
  group_by(name) |> 
  summarise(ssy=median(ssy,na.rm=TRUE),
            cor=cor(value,psi_5, use="pairwise.complete.obs"),
            rmse=exp(sqrt(mean((log(-value)-log(-psi_5))^2,na.rm=TRUE))))

t_test=as.data.frame(psi,xy=TRUE)|>
  pivot_longer(cols=c("psi_2","psi_3","psi_4","psi_6","psi_7","psi_8")) |> 
  mutate(ssy=((psi_5-value)/abs(psi_5))) |> 
  drop_na() |> 
  group_by(name) |> 
  summarise(t.test=t.test(value,psi_5)$p.value)
as.data.frame(psi,xy=TRUE)|>
  rename(real=psi_5,
         `LAI=2`=psi_2,
         `LAI=3`=psi_3,
         `LAI=4`=psi_4,
         `LAI=6`=psi_6,
         `LAI=7`=psi_7,
         `LAI=8`=psi_8
  ) |> 
  pivot_longer(cols=c("LAI=2","LAI=3","LAI=4","LAI=6","LAI=7","LAI=8")) |> 
  mutate(dif=cut((real-value)/1000,
                 breaks=c(-Inf, -2, -1, -0.500,-0.25, 0,0.25,0.500,1,2,Inf))) |> 
  ggplot()+
  geom_tile(aes(x=x,y=y,fill=dif))+
  geom_sf(data=europe$geometry, fill=NA)+
  coord_sf()+
  theme_minimal() +
  theme(axis.title=element_blank(),
        axis.text = element_blank(),
        legend.key.size = unit(0.5,"cm"),
        legend.text = element_text(size=10),
        legend.title = element_text(size=11))+
  scale_fill_brewer(palette="RdYlBu")+
  facet_wrap(~name)+
  labs(fill=unname(TeX("$\\Psi_{min}-\\Psi_{new}(MPa)$")))->figS4
### S5
# idem fig 2 with tar_load(safmarg_spbeta_select)


### S7
db.clim |> 
  filter(presence==1) |> 
  filter(psi_eraday_real>(-10000)) |>
  mutate(psi_eraday_real=psi_eraday_real/1000) |> 
  pivot_longer(cols=c("tmin_era","psi_eraday_real","pet","mat","map")) |> 
  mutate(name=case_when(name=="tmin_era"~"$T_{min}$",
                        name=="psi_eraday_real"~"$\\Psi_{min}$",
                        TRUE~name)) |> 
  group_by(name) |> 
  # mutate(value=scale(value)) |> 
  ungroup() |> 
  group_by(species,name,lt50) |> 
  summarise(`5%`=quantile(value,probs=0.05,na.rm=TRUE),
            `95%`=quantile(value,probs=0.95,na.rm=TRUE)) |> 
  pivot_longer(cols=c("5%","95%"),
               names_to = "Quantile",
               values_to = "val") |> 
  ggplot(aes(lt50,val,color=Quantile))+
  geom_point()+
  geom_smooth(method="gam")+
  scale_color_manual(values=c("darkred","lightblue"))+
  theme_minimal()+
  labs(x="LT50")+
  theme(axis.title.y = element_blank(),
        text = element_text(size=11))+
  facet_wrap(~name,
             labeller=as_labeller(TeX,
                                  default = label_parsed),
             scales="free")->figS7

### S8

db.clim |> 
  filter(presence==1) |> 
  filter(psi_eraday_real>(-10000)) |>
  mutate(psi_eraday_real=psi_eraday_real/1000) |> 
  pivot_longer(cols=c("tmin_era","psi_eraday_real","pet","mat","map")) |> 
  mutate(name=case_when(name=="tmin_era"~"$T_{min}$",
                        name=="psi_eraday_real"~"$\\Psi_{min}$",
                        TRUE~name)) |> 
  group_by(name) |> 
  # mutate(value=scale(value)) |> 
  ungroup() |> 
  group_by(species,name,p50) |> 
  summarise(`5%`=quantile(value,probs=0.05,na.rm=TRUE),
            `95%`=quantile(value,probs=0.95,na.rm=TRUE)) |> 
  pivot_longer(cols=c("5%","95%"),
               names_to = "Quantile",
               values_to = "val") |> 
  ggplot(aes(p50,val,color=Quantile))+
  geom_point()+
  geom_smooth(method="gam")+
  scale_color_manual(values=c("darkred","lightblue"))+
  theme_minimal()+
  labs(x="P50")+
  theme(axis.title.y = element_blank(),
        text = element_text(size=11))+
  facet_wrap(~name,
             labeller=as_labeller(TeX,
                                  default = label_parsed),
             scales="free")->figS8

## S9

predictor=c("overlap_drought","overlap_shade") 
df.mod.overlap=safmarg_sp_select |> 
  left_join(df.species) |> 
  left_join(margin.limit) |> 
  select(taxa,species,
         matches(c("inflex","^t_")),
         hsm.valid.3,fsm.valid.2,
         overlap_drought,overlap_shade,lat.mean,
         mod) |> 
  mutate(across(matches("inflex"),
                ~100-.)) |> 
  pivot_longer(cols=c("overlap_drought","overlap_shade"), 
               names_to = "pred",
               values_to = "pred_val") |> 
  pivot_longer(cols=matches(c("inflex","t_")),
               names_to = "par",
               values_to = "par_val")

df.lm.overlap<-data.frame(pred=rep(predictor,2)) |> 
  crossing(par=c("inflex_fsm","inflex_hsm","t_hsm","t_fsm"),
           mean=NA,
           sd=NA,
           pval=NA,
           conf.up=NA,
           conf.low=NA)


for(preds in predictor){
  print(preds)
  for(pars in c("inflex_fsm","inflex_hsm","t_hsm","t_fsm")){
    print(pars)
    df.mod.pred<-df.mod.overlap |>
      filter(pred==preds) |>
      filter(par==pars) |> 
      filter(!is.na(par_val)&
               !is.na(pred_val)) #|> 
    # mutate(par_val=scale(par_val),
    #        pred_val=scale(pred_val))
    if(grepl("hsm",pars)){
      df.mod.pred<-df.mod.pred |>
        # filter(hsm.valid.3==TRUE) |> 
        filter(mod%in%c("2var","hsm"))|>
        mutate(par_val=scale(par_val),
               pred_val=scale(pred_val))
    }else{
      df.mod.pred<-df.mod.pred |>
        # filter(fsm.valid.2==TRUE) |> 
        filter(mod%in%c("2var","fsm")) |>
        mutate(par_val=scale(par_val),
               pred_val=scale(pred_val))
    }
    summary=as.data.frame(summary(lm(par_val~pred_val,df.mod.pred))$coefficients)
    conf=confint(lm(formula = par_val~pred_val, data = df.mod.pred))
    print(summary(lm(par_val~pred_val,df.mod.pred)))
    df.lm.overlap[df.lm.overlap$pred==preds&
                    df.lm.overlap$par==pars,
                  c("mean","sd","pval","conf.up","conf.low")]=
      c(summary["pred_val",c("Estimate","Std. Error","Pr(>|t|)")],
        conf["pred_val",1],
        conf["pred_val",2])
  }
}


df.lm.overlap<-df.lm.overlap |> 
  mutate(signif1=pval<0.001,
         signif2=pval>0.001&pval<0.01,
         signif3=pval>0.01&pval<0.05)

df.lm.overlap |> 
  # mutate(quant=recode(quant, quant05 = "5%", quant95 = "95%")) |> 
  mutate(pred=case_when(pred=="overlap_drought"~"Competitive dominance  \nfor drought tolerance",
                        pred=="overlap_shade"~"Competitive dominance  \nfor shade tolerance")
  ) |> 
  separate(col = par,
           sep="_", into=c("par","mod")) |> 
  mutate(par=case_when(par=="inflex"~"Inflexion index \n I",
                       par=="t"~"Threshold \n t")) |> 
  ggplot(aes(x=mean,y=par,color=mod))+
  scale_color_manual(values=c("lightblue","darkred"))+
  geom_point(size=2.5,alpha=0.6,position=position_dodge(width=0.5))+ 
  geom_pointrange(aes(xmin=conf.low,xmax=conf.up),position=position_dodge(width=0.5))+ 
  geom_text(aes(label = ifelse(signif3, "*",
                               ifelse(signif2,"**",
                                      ifelse(signif3,"***",""))),
                x=mean),
            size=8,
            show.legend = FALSE,
            vjust=0.75)+
  geom_vline(xintercept = 0,linetype="dotted")+
  theme_minimal()+
  facet_wrap(~pred)+
  theme(axis.title= element_blank(),
        axis.ticks.x = element_line(linewidth = 0.6),
        axis.line.x = element_line(linewidth = 0.6),
        panel.grid = element_blank(),
        # legend.position = "top",
        legend.title=element_blank(),
        text=element_text(size=11),
        axis.text.y = element_text(hjust=0.5)) -> figs9
