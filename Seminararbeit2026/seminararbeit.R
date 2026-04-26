# >NB: -------------------------------------------------------------------------
# - Please make sure to install our package uniLM as described in 
#   https://github.com/Laefritzia/uniLM/README.md
# - documentation to the uniLM functions can be accessed via ?function, 
#    or in the man/ folder in the github.
# 
# - seminararbeit.R contains code to reproduce the plots in our Seminar paper.
# - Simulation results were computed as sees in simulation.R
# 
# Sarah Dvorak, Laetitia Fritz, 25.04.2026
# ------------------------------------------------------------------------------

library(uniLM) # how to install: see https://github.com/Laefritzia/uniLM
library(ggplot2)
library(dplyr)
library(tidyr)
library(flextable)

set.seed(1807)
theme_set(theme_minimal()+theme(legend.position = "bottom"))

# choose scenario A as described in the paper
dataA<-corrupt_data(n=100, scenario="a",O=0.1)

## > Block ---------------------------------------------------------------------

### >> choice of K -------------------------------------------------------------

# this is almost always right, but bc of stochastic blocks may overestimate slightly
#bestK<-adaptK(dataA, algorithm="ADMM",maxiter=20)

true_beta<-dataA$beta
K_grid <- 2:(nrow(dataA$data)/2)

# calculate models for different Ks
b_s <- sapply(K_grid, function(k){
  res<-LM(data=dataA, MOMalgorithm="ADMM", K=k, nboot=0, parallel=FALSE)
  res$MOM$b
})

# euclidean distance between true_beta and estimators
euclid<-apply(b_s, MARGIN=2, function(b){
  sum((b-true_beta)^2) # same rule as in adaptK
})

# choose best one
bestK<-K_grid[which.min(euclid)]
n_outlier<- sum(dataA$data$clean %in% 0)# compare to n outliers

# plot all
K_all<-data.frame(K=K_grid, euclid=euclid) |> 
  ggplot(aes(K, euclid))+
  geom_point(alpha=0.5, color="Darkblue") +
  geom_vline(xintercept=bestK,color="Darkred") +
  annotate("text", 
           x = bestK, 
           y = max(euclid, na.rm = TRUE)/2, 
           label = "bestK", 
           vjust = -0.5, 
           hjust = -0.1, 
           size = 3.5, 
           color = "Darkred")+
  labs(y="Quadratic Error",
       title="Different choices of K", subtitle="Quadratic Error of Estimator")

# plot zoomed
K_zoom <- K_all +
  ylim(c(0,0.1))+
  annotate("text", 
           x = bestK, 
           y = 0.1/2, 
           label = "bestK", 
           vjust = -0.5, 
           hjust = -0.1, 
           size = 3.5, 
           color = "Darkred")

# K_all
# K_zoom


### >> Outlier detection -------------------------------------------------------

mod_a0<-LM(data=dataA, MOMalgorithm="ADMM", nboot=0)

#plot=FALSE stores data
mom_plot<-plot(mod_a0$MOM, plot=FALSE) # object of class MOM: dispatches plot.MOM (mom results)
LM_plot<-plot(mod_a0, plot=FALSE)      # object of class LM: dispatches plot.LM (general diagnostics)

# Outlier detection plots:
plot(mod_a0$MOM, which=3) # mom-block-frequency
plot(mod_a0$MOM, block_p=1/mod_a0$MOM$K, which=3) # mom-block-frequency with custom threshold
plot(mod_a0, which=3) # cooks distance

mom_outs<-which(mom_plot$outlier_data$outlier %in% "Outlier")
cook_outs<-which(LM_plot$cook_d > 2*length(mod_a0$beta)/nrow(mod_a0$data))
true_outs<-which(mod_a0$data$clean == 0)

#mom_outs
#cook_outs
#true_outs
#mom_outs_threshold

#bp<-(1/mod_a0$MOM$K) - sqrt((1-1/mod_a0$MOM$K)*(1/mod_a0$MOM$K)/mod_a0$MOM$n) #minus expected SE (variance of a frequncy: binomial)

### >> Foxed vs. Stochastic blocks  -------------------------------------------

mod_ADMM_r<- LM(dataA, MOMalgorithm=c("ADMM"), stochastic=FALSE, nboot=0)
mod_ADMM_s<- LM(dataA, MOMalgorithm=c("ADMM"), stochastic=TRUE, nboot=0)

# iteration is more wiggly for stochastic
plot(mod_ADMM_r$MOM,which=1:2)
plot(mod_ADMM_s$MOM,which=1:2)

# outlier detection makes no sense with robust blocks
plot(mod_ADMM_r$MOM,which=3)
plot(mod_ADMM_s$MOM,which=3)


## > Simulation Results --------------------------------------------------------

# this is with K outside of Boot
# for reproducitibility consult github.com/Laefritzia/uniLM/Seminararbeit2026/simulation.R
load("Seminar/res100nboot100_H1_Kout_2026-04-26.RData")
load("Seminar/res100nboot100_H0_Kout_2026-04-26.RData")

# 3 Schaetzer, 100 R, 24 Szenarien, laenge2 beta
est<-c("MOM", "OLS", "MEST")

params <- expand.grid(
  scenario=letters[1:4],
  O = c(0.1, 0.3, 0.5),#c(0.1, 0.3, 0.45, 0.5, 0.55, 0.6)
  OO =c(100), #10
  alg = c("GD", "ADMM"),
  nboot = c(100),#500 #1000
  stringsAsFactors = FALSE
)

alpha <- 0.05
Q <- stats::qnorm(1-alpha/2)
true_beta <- c(0,0) #H1 c(1,3), #H0 c(0,0)

bhat<-lapply(est,
  function(e){
    res<-sapply(1:nrow(params), function(x){colMeans(res_b[[e]][,x,])})
    rownames(res) <- c("b0", "b1")
    colnames(res) <- paste0("s",1:nrow(params))
    res
  }) |> 
  setNames(est)

sehat<-lapply(est,
             function(e){
               res<-sapply(1:nrow(params), function(x){colMeans(res_se[[e]][,x,])})
               rownames(res) <- c("se0", "se1")
               colnames(res) <- paste0("s",1:nrow(params))
               res
             }) |> 
  setNames(est)

z <- lapply(est,
            function(e){
              res<-lapply(1:nrow(params), function(x){
                res_b[[e]][,x,] / res_se[[e]][,x,]
                })
              res
            }) |> 
  setNames(est)

bias_mean <- lapply(est,
               function(e){
                 res<-sapply(1:nrow(params), function(x){
                   colMeans(abs(sweep(res_b[[e]][,x,],2,true_beta)))
                 })
                 rownames(res) <- c("bias0", "bias1")
                 colnames(res) <- paste0("s",1:nrow(params))
                 res
               }) |> 
  setNames(est)

bias_median <- lapply(est,
                    function(e){
                      res<-sapply(1:nrow(params), function(x){
                        apply(abs(sweep(res_b[[e]][,x,],2,true_beta)), 2, median)
                      })
                      rownames(res) <- c("bias0", "bias1")
                      colnames(res) <- paste0("s",1:nrow(params))
                      res
                    }) |> 
  setNames(est)

# under H1
pwr<-lapply(est,
              function(e){
                    res<-sapply(1:nrow(params), function(x){
                      colMeans(abs(z[[e]][[x]]) > Q)
                    })
                    rownames(res) <- c("pwr0", "pwr1")
                    colnames(res) <- paste0("s",1:nrow(params))
                    res
                }) |> 
  setNames(est)

# under H0
T1<-lapply(est,
            function(e){
              res<-sapply(1:nrow(params), function(x){
                colMeans(
                  2*(1-stats::pnorm(abs(z[[e]][[x]]))) < alpha
                  )
              })
              rownames(res) <- c("p0", "p1")
              colnames(res) <- paste0("s",1:nrow(params))
              res
            }) |> 
  setNames(est)



all(res_conv$MOM==100)#nboot was 100

# Format for Latex or aggregation
bigtbl<-function(res, caption,b1=TRUE,flextable=TRUE){
  
  ind<-1+as.numeric(b1)#b0 or b1
out<-params |> 
  mutate(
    MOM_b1 =formatC(res$MOM[ind,],digits = 4),
    OLS_b1 =formatC(res$OLS[ind,],digits = 4),
    MEST_b1=formatC(res$MEST[ind,],digits = 4)
  ) |> 
  pivot_wider(
    names_from = alg, 
    values_from = c(MOM_b1, OLS_b1, MEST_b1),
    names_glue = "{alg}_{.value}"
  ) |> 
  select(-c(OO,nboot, matches("GD_OLS|GD_MEST"))) |> 
  rename(OLS_b1=ADMM_OLS_b1, MEST_b1=ADMM_MEST_b1) |> 
  arrange(scenario)

# the Outlier info is meaningless in scenario D, technically we also reran
# the model three times, the fluctuation is just due to epsilon~N(0,1*0.3^2) noise.
d_average<-out |> 
  mutate(across(-scenario, as.numeric))|>
  group_by(scenario) |> 
  mutate(across(-c(O),function(x){
    mean(x)
  })) |>
  ungroup() |>
  mutate(O=0,
         across(-c(scenario,O),~{formatC(.x,digits=1)})) |> 
  slice_tail()

out<-out |> 
  filter(!(scenario%in% "d")) |> 
  bind_rows(d_average)
  
  
if(!b1) {
  out<-out |> rename(OLS_b0=OLS_b1,MEST_b0=MEST_b1,ADMM_MOM_b0=ADMM_MOM_b1,
                     GD_MOM_b0=GD_MOM_b1)
}
if(!flextable) return(out)
  out |> 
  flextable() |>
  merge_v(j=~scenario) |>
  hline(j=1:6,i=c(3,6,9),part="body") |>
  vline(j=c(2),part="body") |>
  align(align = "center", part = "all") |>
  set_caption(caption=caption) |> 
  autofit()
}


b1_bias_mean_H1 <-bigtbl(bias_mean, "Simulation Results: Bias (H1).")
b1_bias_median_H1<-bigtbl(bias_median, "Simulation Results: Bias (H1) - Median of Simulations")
b1_pwr <-bigtbl(pwr, "Power (alpha=0.05)")
b1_T1 <- bigtbl(T1, "Type I Error (alpha=0.05)")

b0_bias_mean_H1<-bigtbl(bias_mean, "Simulation Results: Bias (H1)", b1=FALSE)
b0_bias_median_H1<-bigtbl(bias_median, "Simulation Results: Bias (H1)", b1=FALSE)
b0_pwr<-bigtbl(pwr, "Power (alpha=0.05)", b1=FALSE)
b0_T1<-bigtbl(T1, "Type I Error (alpha=0.05)",b1=FALSE)

# slope pwr scenario a 
bigtbl(pwr, "",flextable=FALSE) |> 
  filter(scenario %in% "a") |> 
  mutate(across(-scenario, as.numeric))|> 
  pivot_longer(cols=-c(scenario,O)) |> 
  mutate(estimate=factor(name)) |> 
  ggplot(aes(O, value,color=estimate)) + 
  scale_color_manual(values = c("#009E73","#CC79A7","#E69F00","#56B4E9"))+
  geom_line()+geom_point()+facet_wrap(~scenario)+
  labs(title="Power ~ Outlier Percentage O",subtitle="100 Simulations, Scenario a.",
       y="pwr")

# slope mean bias scenario a without GD
bigtbl(bias_mean, "", flextable = FALSE) |> 
  select(-GD_MOM_b1) |> 
  filter(scenario %in% "a") |> 
  mutate(across(-scenario, as.numeric)) |> 
  pivot_longer(cols=-c(scenario,O)) |> 
  mutate(estimate=factor(name)) |> 
  ggplot(aes(O,value,color=estimate,group=estimate)) +
  scale_color_manual(values = c("#009E73","#E69F00","#56B4E9"))+
  geom_line() +geom_point()+
  labs(title="Mean Bias ~ Outlier Percentage O",subtitle="100 Simulations, Scenarios a",
       y="Bias")

#slope pwr
# bigtbl(pwr, "",flextable=FALSE) |> 
#   mutate(across(-scenario, as.numeric))|> 
#   pivot_longer(cols=-c(scenario,O)) |> 
#   mutate(estimate=factor(name)) |> 
#   ggplot(aes(O, value,color=estimate)) + 
#   geom_line()+geom_point()+facet_wrap(~scenario)+
#   labs(title="Power ~ Outlier Percentage O",subtitle="100 Simulations, Scenarios a-d.",
#        y="pwr")

#> Appendix: testing reshuffling of blocks ONCE/TWICE------------------------

# Lerasle and Lecué (2020) propose to shuffle random blocks TWICE per iterations,
# so once for each of the two candidate g and f. We implemented only ONE resampling
# at the start of each iteration (so the same assimgment for g and f).
# This is computationally faster.

mods_shuffle1<-list()
for(i in 1:100){
  res<-LM(dataA, MOMalgorithm = "ADMM",stochastic=TRUE,nboot=0,parallel=FALSE)
  mods_shuffle1[[i]]<-res$MOM
}

# here we duplicated the the shuffling logic inside 'uniLM:::calculate_mom'
# and loaded the temporary manipulated function into our R-session.
mods_shuffle2<-list()
for(i in 1:100){
  res<-LM(dataA, MOMalgorithm = "ADMM",stochastic=TRUE,nboot=0,parallel=FALSE)
  mods_shuffle2[[i]]<-res$MOM
}

plot_its<-function(mod, var){
purrr::map_dfr(mod, function(x) tibble::tibble(y=x[[var]],i=1:100)) |> 
  mutate(run=factor(rep(1:100,each=100)))|> 
  group_by(i) |> 
  mutate(med=median(y)) |> 
  ungroup() |> 
  ggplot(aes(i, y, group=run))+
  geom_line(color="Darkblue")+
  geom_line(aes(i, med,color="Median"),lwd=1.2)+
  labs(y=var,x="iterations",color="")
}

plot_its(mods_shuffle1, "mom_obj") +
  labs(title="MOM-objective trajectory of 100 regressions", subtitle="Block resampling ONCE per iteration")
plot_its(mods_shuffle1, "mom_err")+
  labs(title="Estimation Error trajectory of 100 regressions", subtitle="Block resampling ONCE per iteration")
plot_its(mods_shuffle2, "mom_obj") +
  labs(title="MOM-objective trajectory of 100 regressions", subtitle="Block resampling TWICE per iteration")
plot_its(mods_shuffle2, "mom_err") +
  labs(title="Estimation Error trajectory of 100 regressions", subtitle="Block resampling TWICE per iteration")

