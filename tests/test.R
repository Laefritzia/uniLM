# > testing MOM  -----------------
if(FALSE){

  temp<-corrupt_data(n=100, O=0.1) # bis zu 0.2 gehts noch ganz gut, mit 0.1 besser
  data<-temp$data; beta<-temp$beta; rm(temp)

  # number of blocks K:
  # lepski K: the paper mentions that K (number of blocks)
  # has to be at least 8 bigger than the number of outliers (but
  # this number can be unknown) to be optimal. they present an adaptive choice
  # for K, I didn't really understand it yet, so im looking through all

  K_grid<-2:(nrow(data)/2)

  # >> Subgradient -----------------
  res<- sapply( # looking for best K
    lapply(K_grid, function(K){
      MOM.LM(data, beta, K, algorithm="GD",maxiter=20)
    }),
    function(x){
      sum((x$b-beta)^2)
    })
  #sum((data$clean%in%0)*8);sum((data$clean%in%0)*2)
  #K<-K_grid[which.min(res)]
  r1<-MOM.LM(data, beta, K_grid[which.min(res)], algorithm="GD",maxiter=100)
  plot(0:100, r1$mom_obj)
  plot(0:100, r1$mom_err)
  abs(r1$b-beta)

  clean<-data$clean %in% 1
  ols<-with(data,{
    X<-model.matrix(as.formula("y~x2"),data=data)
    solve(t(X)%*%X)%*%t(X)%*%y
  })
  ols_clean<-lm("y~x2",data=data[data$clean==1,])
  r1$b; ols; ols_clean$coefficients
  plot(data$x2, data$y)
  abline(a=r1$b[1], b=r1$b[2],col="green") # unser MOM
  abline(a=ols_clean[1], b=ols_clean[2],col="yellow") # clean OLS (~truth)
  abline(a=ols[1], b=ols[2],col="red") # naiver OLS (schlecht)
  plot(data$x2[clean], data$y[clean])
  abline(a=r1$b[1], b=r1$b[2],col="green")
  abline(a=ols_clean[1], b=ols_clean[2],col="yellow")
  abline(a=ols[1], b=ols[2],col="red") #lol

  # >> ADMM -----------------
  res2<- sapply( # looking for best K
    lapply(K_grid, function(K){
      MOM.LM(data, beta, K, algorithm="ADMM",maxiter=20)
    }),
    function(x){
      sum((x$b-beta)^2)
    })
  r2<-MOM.LM(data, beta, K_grid[which.min(res2)], algorithm="ADMM",maxiter=100)
  plot(0:100, r2$mom_obj)
  plot(0:100, r2$mom_err)
  abs(r2$b-beta)


  r2$b; ols; ols_clean$coefficients
  plot(data$x2, data$y)
  abline(a=r2$b[1], b=r2$b[2],col="green") # unser MOM
  abline(a=ols_clean[1], b=ols_clean[2],col="yellow") # clean OLS (~truth)
  abline(a=ols[1], b=ols[2],col="red") # naiver OLS (schlecht)
  plot(data$x2[clean], data$y[clean])
  abline(a=r2$b[1], b=r2$b[2],col="green")
  abline(a=ols_clean[1], b=ols_clean[2],col="yellow")
  abline(a=ols[1], b=ols[2],col="red") #lol


}
