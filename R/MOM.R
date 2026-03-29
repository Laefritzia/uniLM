
# > subgradient descent  -----------------
subgrad<-function(b,data,form=as.formula("y~x2")){
  y<-data[,all.vars(form)[1]]
  X<-model.matrix(form,data=data)
  t(X)%*%(X%*%b-y)
}

# > admm descent/ascent  -----------------

# weil kein lasso brauchen wir soft treshhold nicht
admm<-function(rho,z,u,data,form=as.formula("y~x2")){
  y<-data[,all.vars(form)[1]]
  X<-model.matrix(form,data=data)

  b<-solve(t(X)%*%X+rho*diag(ncol(X)), #kxk
           t(X)%*%y+rho*z-u) #kx1 *1xk
  z<-b+ u/rho
  u<-u+rho*(b-z)
  list(b=b, z=z, u=u)
}

# > median block  -----------------
med_block <- function(data, form=as.formula("y~x2"),
                       K, b, b_prime){

  y<-data$y
  X<-model.matrix(form,data=data)
  blocks<-sample(factor(rep(1:K, length.out=nrow(data))))

  # for each block we calculate mean loss
  means_loss<-sapply(1:K, function(B){

    Xk <- X[blocks%in%B,,drop=FALSE]
    yk <- y[blocks%in%B]

    sum((Xk%*%b-yk)^2) - sum((Xk%*%b_prime-yk)^2)

  })

  # choose block which is closest to median
  med_ind<-which.min(abs(median(means_loss)-means_loss))[1]

  return(data[blocks %in%med_ind,])

}

# > MOM.LM -----------------
MOM.LM <- function(data, beta, K,
                   form=as.formula("y~x2"),
                   algorithm=c("GD", "ADMM"),
                   stepsize=0.01,
                   maxiter=100,
                   tol=10^-6){

data<-data[,all.vars(form)]
X<-model.matrix(form,data=data)
y<-data[,all.vars(form)[1]]

alg<-match.arg(algorithm)

# > Gradient Descent -----------------
if(alg%in% "GD"){ # Gradient Descent

b<-b_prime<-rep(0,ncol(X))
mom_obj<-mom_err<-numeric(maxiter)

iter<-0
while(TRUE) {
iter<-iter+1
 # medium worst block:(maximization)
  block<- med_block(data,form,K,b,b_prime)
 # gradient descent:(minimization)
  b    <- b - (stepsize)*subgrad(b,form,data=block)#/sqrt(iter)
 # same with new b for b_prime
  block<- med_block(data,form,K,b,b_prime)
  b_prime<- b_prime - (stepsize)*subgrad(b_prime,form,data=block)#/sqrt(iter)

  #mom(l_b - l_b_prime)
  mom_obj[iter]<- sum((X%*%b-y)^2) - sum((X%*%b_prime-y)^2)
  mom_err[iter]<- sqrt(sum((b-beta)^2))

  if(iter>maxiter) break
}

# maybe give MOM class for plots?
return(list(
  b=b,
  mom_obj=mom_obj,
  mom_err=mom_err
))
}

# > ADMM -----------------
if(alg%in% "ADMM"){

 b<-b_prime<-u<-u_prime<-z<-z_prime<-rep(0,ncol(X))
 mom_obj<-mom_err<-numeric(maxiter)
 rho<-5
 iter<-0
 while(TRUE) {
   iter<-iter+1
   # DESCENT
   block<- med_block(data,form,K,b,b_prime)
   admmD <- admm(rho,z,u,data=block)
   b<-admmD$b; z<- admmD$z; u<- admmD$u
   # ASCENT
   block<- med_block(data,form,K,b,b_prime)
   admmA <- admm(rho,z_prime,u_prime,data=block)
   b_prime<-admmA$b; z_prime<- admmA$z; u_prime<- admmA$u

   #mom(l_b - l_b_prime)
   mom_obj[iter]<- sum((X%*%b-y)^2) - sum((X%*%b_prime-y)^2)
   mom_err[iter]<- sqrt(sum((b-beta)^2))

   if(iter>maxiter) break
 }

 # maybe give MOM class for plots?
 return(list(
   b=b,
   mom_obj=mom_obj,
   mom_err=mom_err
 ))

}

}




