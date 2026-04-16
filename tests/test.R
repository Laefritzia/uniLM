
# > testing MOM  -----------------
if(FALSE){
  devtools::load_all()

  data<-corrupt_data(n=100, scenario="a",O=0.1) # bis zu 0.2 gehts noch ganz gut, mit 0.1 besser

  # number of blocks K:
  # lepski K: the paper mentions that K (number of blocks)
  # has to be at least 8 bigger than the number of outliers (but
  # this number can be unknown) to be optimal. they present an adaptive choice
  # for K, I didn't really understand it yet, so im looking through all

  K_grid<-2:(nrow(data$data)/2)

  # >> Subgradient -----------------
  res<- sapply( # looking for best K
    lapply(K_grid, function(K){
      MOM(data, K=K, algorithm="GD",maxiter=20,stochastic=FALSE)
    }),
    function(x){
      sum((x$b-data$beta)^2)
    })
  res_s<- sapply( # looking for best K
    lapply(K_grid, function(K){
      MOM(data, K=K, algorithm="GD",maxiter=20,stochastic=TRUE)
    }),
    function(x){
      sum((x$b-data$beta)^2)
    })
  #sum((data$data$clean%in%0)*8);sum((data$data$clean%in%0)*2)
  #K<-K_grid[which.min(res)]
  r1<-MOM(data, K=K_grid[which.min(res)], algorithm="GD",maxiter=100,stochastic=FALSE)
  r1_s<-MOM(data, K=K_grid[which.min(res_s)], algorithm="GD",maxiter=100,stochastic=TRUE)
  plot(1:100, r1$mom_obj)
  plot(1:100, r1$mom_err)
  plot(1:100, r1_s$mom_obj)
  plot(1:100, r1_s$mom_err)
  abs(r1$b-data$beta)
  abs(r1_s$b-data$beta)
  clean<-data$data$clean %in% 1

  boxplot(scores~clean,data.frame(scores=r1$scores, clean=clean))
  boxplot(scores~clean,data.frame(scores=r1_s$scores, clean=clean))
  mean(r1$scores[clean]); mean(r1$scores[!clean])
  mean(r1_s$scores[clean]); mean(r1_s$scores[!clean])

  sum(r1$scores[clean]== 0)
  sum(r1_s$scores[clean]== 0)

  M<-matrix(0,nrow=10,ncol=6)
  colnames(M)<-c(apply(
    expand.grid(c("rob", "stoch"), c("_CleanScore", "_UncleanScore")),
    1,paste,collapse=""),
    "rob_CleanIgnore", "stoch_CleanIgnore")
  for(i in 1:nrow(M)){
  r1<-MOM(data, K=K_grid[which.min(res)], algorithm="GD",maxiter=100,stochastic=FALSE)
  r1_s<-MOM(data, K=K_grid[which.min(res_s)], algorithm="GD",maxiter=100,stochastic=TRUE)

  M[i,]<-cbind(
    mean(r1$scores[clean]), mean(r1_s$scores[clean]),
    mean(r1$scores[!clean]), mean(r1_s$scores[!clean]),
    round(sum(r1$scores[clean]== 0)), round(sum(r1_s$scores[clean]== 0))
    )
  }
  colMeans(M)

  ols<-with(data$data,{
    X<-model.matrix(as.formula("y~x2"),data=data$data)
    solve(t(X)%*%X)%*%t(X)%*%y
  })
  ols_clean<-lm("y~x2",data=data$data[data$data$clean==1,])
  r1$b; r1_s$b; ols; ols_clean$coefficients
  plot(data$data$x2, data$data$y)
  abline(a=r1$b[1], b=r1$b[2],col="green") # unser MOM (non stochastic)
  abline(a=r1_s$b[1], b=r1_s$b[2],col="darkgreen") # unser MOM (stochastic)
  abline(a=ols_clean[1], b=ols_clean[2],col="yellow") # clean OLS (~truth)
  abline(a=ols[1], b=ols[2],col="red") # naiver OLS (schlecht)
  plot(data$data$x2[clean], data$data$y[clean])
  abline(a=r1$b[1], b=r1$b[2],col="green")
  abline(a=r1_s$b[1], b=r1_s$b[2],col="darkgreen") # unser MOM (stochastic)
  abline(a=ols_clean[1], b=ols_clean[2],col="yellow")
  abline(a=ols[1], b=ols[2],col="red") #lol

  # >> ADMM -----------------
  res2<- sapply( # looking for best K
    lapply(K_grid, function(K){
      MOM(data, K=K, algorithm="ADMM",maxiter=20,stochastic=FALSE)
    }),
    function(x){
      sum((x$b-data$beta)^2)
    })
  res2_s<- sapply( # looking for best K
    lapply(K_grid, function(K){
      MOM(data, K=K, algorithm="ADMM",maxiter=20,stochastic=TRUE)
    }),
    function(x){
      sum((x$b-data$beta)^2)
    })
  r2<-MOM(data, K=K_grid[which.min(res2)], algorithm="ADMM",maxiter=100,stochastic=FALSE)
  r2_s<-MOM(data, K=K_grid[which.min(res2_s)], algorithm="ADMM",maxiter=100,stochastic=TRUE)
  plot(1:100, r2$mom_obj)
  plot(1:100, r2$mom_err)
  abs(r2$b-data$beta)
  plot(1:100, r2_s$mom_obj)
  plot(1:100, r2_s$mom_err)
  abs(r2_s$b-data$beta)

  boxplot(scores~clean,data.frame(scores=r2$scores, clean=clean))
  boxplot(scores~clean,data.frame(scores=r2_s$scores, clean=clean))
  mean(r2$scores[clean]); mean(r2$scores[!clean])
  mean(r2_s$scores[clean]); mean(r2_s$scores[!clean])

  sum(r2$scores[clean]== 0)
  sum(r2_s$scores[clean]== 0)

  M<-matrix(0,nrow=10,ncol=6)
  colnames(M)<-c(apply(
    expand.grid(c("rob", "stoch"), c("_CleanScore", "_UncleanScore")),
    1,paste,collapse=""),
    "rob_CleanIgnore", "stoch_CleanIgnore")
  for(i in 1:nrow(M)){
    r2<-MOM(data, K=K_grid[which.min(res2)], algorithm="ADMM",maxiter=100,stochastic=FALSE)
    r2_s<-MOM(data, K=K_grid[which.min(res2_s)], algorithm="ADMM",maxiter=100,stochastic=TRUE)

    M[i,]<-cbind(
      mean(r2$scores[clean]), mean(r2_s$scores[clean]),
      mean(r2$scores[!clean]), mean(r2_s$scores[!clean]),
      round(sum(r2$scores[clean]== 0)), round(sum(r2_s$scores[clean]== 0))
    )
  }
  colMeans(M)



  r2$b; ols; ols_clean$coefficients
  plot(data$data$x2, data$data$y)
  abline(a=r2$b[1], b=r2$b[2],col="green") # unser MOM
  abline(a=r2_s$b[1], b=r2_s$b[2],col="darkgreen") # unser MOM(stochastic)
  abline(a=ols_clean[1], b=ols_clean[2],col="yellow") # clean OLS (~truth)
  abline(a=ols[1], b=ols[2],col="red") # naiver OLS (schlecht)
  plot(data$data$x2[clean], data$data$y[clean])
  abline(a=r2$b[1], b=r2$b[2],col="green")
  abline(a=r2_s$b[1], b=r2_s$b[2],col="darkgreen") # unser MOM(stochastic)
  abline(a=ols_clean[1], b=ols_clean[2],col="yellow")
  abline(a=ols[1], b=ols[2],col="red") #lol


# > plot LM  -----------------

  data<-corrupt_data(n=100, scenario="a",O=0.1) # bis zu 0.2 gehts noch ganz gut, mit 0.1 besser

  plot(data,which=3)
  plot(data) # macht vier plots hintereinander
  plot(data,which=5) # fehler
  plot_data<-plot(data,plot=FALSE) # nur die daten mit plot FALSE

}

