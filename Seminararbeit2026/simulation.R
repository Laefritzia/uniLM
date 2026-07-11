# >NB: -------------------------------------------------------------------------
# - Please make sure to install our package uniLM as described in 
#   https://github.com/Laefritzia/uniLM/README.md
# - documentation to the uniLM functions can be accessed via ?function, 
#    or in the man/ folder in the github.
# 
# Simulation results as seen in main section in Seminarpaper: 
# - adapted K outside of bootstrap.
# - Run separetly for each Hypothesis (mybeta).

# The set-up of this simulation loop and logic in R was heavily supported 
# by the use of a public AI-model (Google: Gemini AI, April 2026).
# 
# Sarah Dvorak, Laetitia Fritz, 25.04.2026
# ------------------------------------------------------------------------------

library(future)
library(uniLM)

set.seed(1807)
R<-100
workers<-parallel::detectCores()-1

plan(sequential)#default
plan(multisession, workers=workers)


# combinations of all the scenarios,params,...etc.
params <- expand.grid(
  scenario=letters[1:4],
  O = c(0.1, 0.3, 0.5),#c(0.1, 0.3, 0.45, 0.5, 0.55, 0.6)
  OO =c(100), #10
  alg = c("GD", "ADMM"),
  nboot = c(100),#500 #1000
  stringsAsFactors = FALSE
)

# algorithm indexes (nboot and algorithm info)
# since data stays the same for the nboo/algorithm combinations,
# we split by the data parameters. Eg. the meaning of the following line:
# the following indexes give the info on running the algorithm 4 times, for the
# data with scenario a, O=0.45, OO=100.
# $a.0.45.100
# [1]  33  81 129 177
dta_alg<-split(1:nrow(params),
      list(params$scenario,params$O, params$OO))

est<-c("MOM", "OLS", "MEST")

# run two times for H0, H1, so its more clear
mybeta<-c(0,0) # c(1,3)

# res_b and res_se are lists respectively, with three arrays for each estimate.
# array has three dimensions: (other way than in python)
# row-col-group, for us: R-param/combination-coefficient(intercept/slope)
res_b<-lapply(est, function(x){
  array(0, c(R, nrow(params), length(mybeta)))
}) |>
  setNames(est)

res_se<-lapply(est, function(x){
  array(0, c(R, nrow(params), length(mybeta)))
}) |>
  setNames(est)

# convergence info
res_conv<-lapply(est, function(x){
  array(0, c(R, nrow(params))) 
}) |>
  setNames(est)

message("Start: ", Sys.time())
for (i in 1:R){
if(i%%5==0) message("Progress: Starting", i, " -th Iteration at ", Sys.time())
  
  for(alg in dta_alg){
    
    dta<-params[alg[1],]
    
    simdata<-corrupt_data(n=100, beta=mybeta,
                          scenario=dta$scenario, O=dta$O, OO=dta$OO)
    
    for (j in alg){
      
      bestK<-adaptK(simdata, algorithm=params$alg[j])
      
      simfit<-LM(data=simdata, K=bestK,
                 MOMalgorithm=params$alg[j],
                 nboot=params$nboot[j],
                 parallel=TRUE)
      
      for (x in est){
        res_b[[x]][i,j,] <-simfit[[x]]$b
        res_se[[x]][i,j,]<-simfit[[x]]$se
      }
      
      res_conv$MOM[i,j]<-simfit$MOM$nboot_conv
    }
    
  }
  
if(i%%10==0) save.image(file="tmp.RData") 
}

message("End: ", Sys.time())
plan(sequential) # to default again
save(res_b, res_se, res_conv,file=paste0("res", R, "nboot100_H0_Kout_", Sys.Date(), ".RData"))

#res_b$MOM[1,,]
#res_b$MOM[,1,]
#res_se$MOM[,1,]
