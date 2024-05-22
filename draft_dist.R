# data.frame(kint=c(rbeta(10000,4,12),rbeta(10000,20,60)),
#            lambda=c(rgamma(10000,6,2),rpareto(10000,5,7)),
#            type=c(rep("a",10000),rep("b",10000))) |>
data.frame(kint=rbeta(10000,1.5,15),
           lambda=rgamma(10000,18,2)) |>
# data.frame(kint=rbeta(10000,20,60),
#            lambda=rpareto(10000,5,7)) |>
  rowwise() |> 
  mutate(ks=rbeta(1,lambda*kint,lambda*(1-kint))) |> 
  # group_by(type) |> summarise_all(.funs = mean)
  summary()
  pivot_longer(cols=everything()) |> 
  ggplot()+
  geom_density(aes(value,color=name))+
  # geom_histogram(data=data.frame(kfix=as.numeric(summary(fit.allsp)$summary[1:38,1])),
  #              aes(kfix),
  #              alpha=0.2)+
  facet_wrap(~name,scales="free")
# K_int ~ beta(2.5,15);
# lambda ~ gamma(6,0.5);
# 
# 
# K_sp ~ beta(lambda * K_int, lambda * (1 - K_int));


data.frame(b_1=rgamma(10000,1.5,0.5),
           b_2=rgamma(10000,6,0.5)) |> 
  rowwise() |> 
  mutate(ks=rbeta(1,b_1,b_2)) |> 
  # summary()
ggplot(aes(ks))+
  geom_density()

b_1 ~ gamma(1.5,0.5);
b_2 ~ gamma(6,0.5);

load("mod_random_binom_fixed_alldata.RData")
plot(density(summary(fit.allsp)$summary[1:38,1]))
