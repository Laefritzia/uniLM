
# > gradient descent  -----------------
gdesc<-function(X,y,b){
  t(X)%*%(X%*%b-y)
}

# > admm descent/ascent  -----------------

# weil kein lasso brauchen wir soft treshhold nicht
admm<-function(X,y,rho,z,u){

  b<-solve(t(X)%*%X+rho*diag(ncol(X)), #kxk
           t(X)%*%y+rho*z-u) #kx1 *1xk
  z<-b+ u/rho
  u<-u+rho*(b-z)
  list(b=b, z=z, u=u)
}

# > median block  -----------------

#' chooses block based on MOM: median of means
#'
#' @param X model.matrix containing covariates and intercept.
#' @param y numeric response vector
#' @param K number of blocks
#' @param b coefficient candidate
#' @param b_prime another coefficient candidate
#' @returns named list with X and y of median block
#' @export
med_block <- function(X,y,K,b,b_prime){

  blocks<-sample(factor(rep(1:K, length.out=nrow(X))))

  # for each block we calculate mean loss
  means_loss<-sapply(1:K, function(B){
    Xk <- X[blocks%in%B,,drop=FALSE]
    yk <- y[blocks%in%B]

    sum((Xk%*%b-yk)^2) - sum((Xk%*%b_prime-yk)^2)
  })

  # choose block which is closest to median
  med_ind<-which.min(abs(stats::median(means_loss)-means_loss))[1]

  return(list(
    X=X[blocks %in%med_ind,,drop=FALSE],
    y=y[blocks %in%med_ind]
  ))

}

# > MOM.LM -----------------

#' robust MOM-estimator for linear regression
#'
#' @param data data.frame of class LM
#' @param K integer specifying number of blocks for MOM-algorithm
#' @param algorithm specifying how to get coefficients, GD (gradient-descent), ADMM(ascent-descent)
#' @returns named list with final b (coefficients), iterative objectives and errors
#' @export
MOM.LM <- function(data, K, algorithm=c("GD", "ADMM"),
                   stepsize=0.01,
                   maxiter=100){

# Checks:
if(!inherits(data, "LM")){
    stop("Data must be of class LM (eg 'uniLM::LM()',uniLM::corrupt_data()')")
}
invisible(validate_LM(data))
invisible(mapply(
  check_1num, list(K, stepsize, maxiter), c("int", "num", "int"))
  )


# parameters:
beta<-data$beta
form<-data$form
data<-data$data[,all.vars(form)]
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
  block<- med_block(X,y,K,b,b_prime)
 # gradient descent:(minimization)
  b    <- b - (stepsize)*gdesc(block$X, block$y, b)#/sqrt(iter)
 # same with new b for b_prime
  block<- med_block(X,y,K,b,b_prime)
  b_prime<- b_prime - (stepsize)*gdesc(block$X, block$y, b_prime)#/sqrt(iter)

  #mom(l_b - l_b_prime)
  mom_obj[iter]<- sum((X%*%b-y)^2) - sum((X%*%b_prime-y)^2)
  mom_err[iter]<- sqrt(sum((b-beta)^2))

  if(iter>maxiter) break
}

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
 rho<-5 # same as in their paper?
 iter<-0
 while(TRUE) {
   iter<-iter+1
   # DESCENT
   block<-  med_block(X,y,K,b,b_prime)
   admmD <- admm(block$X, block$y, rho,z,u)
   b<-admmD$b; z<- admmD$z; u<- admmD$u
   # ASCENT
   block<- med_block(X,y,K,b,b_prime)
   admmA <- admm(block$X, block$y,rho,z_prime,u_prime)
   b_prime<-admmA$b; z_prime<- admmA$z; u_prime<- admmA$u

   #mom(l_b - l_b_prime)
   mom_obj[iter]<- sum((X%*%b-y)^2) - sum((X%*%b_prime-y)^2)
   mom_err[iter]<- sqrt(sum((b-beta)^2))

   if(iter>maxiter) break
 }

 return(list(
   b=b,
   mom_obj=mom_obj,
   mom_err=mom_err
 ))

}

}




