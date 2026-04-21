
# > gradient descent  -----------------
gdesc<-function(X,y,b){
  t(X)%*%(X%*%b-y)
}

# > admm descent/ascent  -----------------
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
#' @param blocks supply block-partition of the data (eg fixed blocks or reshuffled each iteration)
#' @examples
#'  data<-uniLM::corrupt_data(n=100,scenario="a",O=0.1)
#'  X<-model.matrix(data$formula,data=data$data)
#   y<-data$data[,all.vars(data$formula)[1]])
#
#'  med_block(X=X,y=y,K=20,
#'  b=rep(0,length(data$beta)),
#'  b_prime=data$beta
#'  )
#' @returns named list with X and y of median block
#' @export
med_block <- function(X,y,K,b,b_prime, blocks=NULL){

  blocks<-if (is.null(blocks)) sample(factor(rep(1:K, length.out=n))) else blocks

  # for each block we calculate mean loss
  means_loss<-sapply(1:K, function(B){
    Xk <- X[blocks%in%B,,drop=FALSE]
    yk <- y[blocks%in%B]

    sum((Xk%*%b-yk)^2)/nrow(Xk) - sum((Xk%*%b_prime-yk)^2)/nrow(Xk)
  })

  # choose block which is closest to median
  med_ind<-which.min(abs(stats::median(means_loss)-means_loss))[1]

  return(list(
    X=X[blocks %in%med_ind,,drop=FALSE],
    y=y[blocks %in%med_ind],
    med_ind=med_ind
  ))

}


# > Calculate MOM -----------------------------------
calculate_mom <- function(data, formula, beta, K, algorithm=c("GD", "ADMM"),
                          stochastic=TRUE, stepsize=0.01, maxiter=100, seed=NULL){

  # checks and paramaters

  invisible(mapply(
    check_1num, list(K, stepsize, maxiter), c("int", "num", "int"))
  )

  alg<-match.arg(algorithm)
  X<-model.matrix(formula,data=data)
  y<-data[,all.vars(formula)[1]]
  n<-nrow(X)

  if(!is.null(seed)) set.seed(seed)
  fix_blocks<-if(stochastic) NULL else sample(factor(rep(1:K, length.out=n)))

  # > Gradient Descent -----------------
  if(alg%in% "GD"){ # Gradient Descent

    b<-b_prime<-numeric(ncol(X))
    mom_obj<-mom_err<-numeric(maxiter)
    scores<-numeric(n)

    iter<-0
    while(TRUE) {
      iter<-iter+1

      blocks<- if(stochastic) sample(factor(rep(1:K, length.out=n))) else fix_blocks

      # medium worst block:(maximization)
      block<- med_block(X,y,K,b,b_prime, blocks=blocks)
      scores[blocks==block$med_ind] <- scores[blocks==block$med_ind]+1
      # gradient descent:(minimization)
      b    <- b - (stepsize)*gdesc(block$X, block$y, b)#/sqrt(iter)
      # same with new b for b_prime
      block<- med_block(X,y,K,b,b_prime, blocks=blocks)
      scores[blocks==block$med_ind] <- scores[blocks==block$med_ind]+1
      b_prime<- b_prime - (stepsize)*gdesc(block$X, block$y, b_prime)#/sqrt(iter)

      #mom(l_b - l_b_prime)
      mom_obj[iter]<- sum((X%*%b-y)^2)/n - sum((X%*%b_prime-y)^2)/n
      mom_err[iter]<- sqrt(sum((b-beta)^2))

      if(iter>=maxiter) break
    }

    return(list(
      b=b,
      mom_obj=mom_obj,
      mom_err=mom_err,
      scores=scores/maxiter

    ))
  }

  # > ADMM -----------------
  if(alg%in% "ADMM"){

    b<-b_prime<-u<-u_prime<-z<-z_prime<-numeric(ncol(X))
    mom_obj<-mom_err<-numeric(maxiter)
    scores<-numeric(n)
    rho<-5 # same as in their paper?
    iter<-0
    while(TRUE) {
      iter<-iter+1
      blocks<- if(stochastic) sample(factor(rep(1:K, length.out=n))) else fix_blocks
      # DESCENT
      block<-  med_block(X,y,K,b,b_prime,blocks=blocks)
      scores[blocks==block$med_ind] <- scores[blocks==block$med_ind]+1
      admmD <- admm(block$X, block$y, rho,z,u)
      b<-admmD$b; z<- admmD$z; u<- admmD$u
      # ASCENT
      block<- med_block(X,y,K,b,b_prime,blocks=blocks)
      scores[blocks==block$med_ind] <- scores[blocks==block$med_ind]+1
      admmA <- admm(block$X, block$y,rho,z_prime,u_prime)
      b_prime<-admmA$b; z_prime<- admmA$z; u_prime<- admmA$u

      #mom(l_b - l_b_prime)
      mom_obj[iter]<- sum((X%*%b-y)^2)/n - sum((X%*%b_prime-y)^2)/n
      mom_err[iter]<- sqrt(sum((b-beta)^2))

      if(iter>=maxiter) break
    }

    return(list(
      b=b,
      mom_obj=mom_obj,
      mom_err=mom_err,
      scores=scores/maxiter

    ))
  }

}

# > MOM ----------------

#' robust MOM-estimator for linear regression
#'
#' @param data data.frame containing regression data
#' @param formula formula specifiyng the relationship of response and predictors in data
#' @param beta true coefficient for relationship (for simulation purposes)
#' @param K integer specifying number of blocks for MOM-algorithm
#' @param algorithm specifying how to get coefficients, GD (gradient-descent), ADMM(ascent-descent)
#' @param stochastic TRUE/FALSE decides, if blocks should be shuffled randomly in every iteration (stochastic==TRUE)
#'                   or based on a fixed starting partition (stochastic==FALSE)
#' @param nboot for variance estimation: number of bootstrap iterations (with resampling blocks...stochastic==TRUE)
#' @param ... optional arguments like stepsize, maxiter, seed, ...
#' @returns named list with final b (coefficients), iterative objectives and errors
#' @examples
#'   data<-uniLM::corrupt_data(n=100,scenario="a",O=0.1)
#'   res<-MOM(data, K=12, algorithm="GD")
#'   plot(0:100, res$mom_obj)
#'   plot(0:100, res$mom_err)
#'   abs(res$b-true_beta)
#' @export
MOM <- function(data, formula=NULL, beta=NULL, K, algorithm=c("GD", "ADMM"),
                stochastic=TRUE, nboot=0,...){
  UseMethod("MOM")
}


#' @rdname MOM
#' @export
MOM.default <- function(data, formula, beta, K, algorithm=c("GD", "ADMM"),
                        stochastic=TRUE, nboot=0, ...){

  # checks:

  invisible(validate_LM(new_LM(list(data=data,formula=formula,beta=beta))))

  # bootstrap logic:

  boots<-list(); length(boots) <- nboot

  if (nboot>0){
    for (i in 1:nboot){
      boots[[i]]<-calculate_mom(data=data, formula=formula, beta=beta, K=K,
                                algorithm=algorithm, stochastic=stochastic, ...)
    }

    return(list(
      betas=sapply(boots, function(x) x$b),
      b=rowMeans(betas),
      vcov = cov(t(betas)),
      se=sqrt(diag(vcov)),
      scores=rowMeans(sapply(boots,function(x)x$scores))
    ))

  } else{

    calculate_mom(data=data, formula=formula, beta=beta, K=K,
                  algorithm=algorithm, stochastic=stochastic, ...)
  }
}

#' @rdname MOM
#' @export
MOM.LM <- function(data, formula=NULL, beta=NULL, K, algorithm=c("GD", "ADMM"),
                   stochastic=TRUE, nboot=0,...){

  # checks and structure:

  if(!inherits(data, "LM")){
    stop("Data must be of class LM (eg 'uniLM::LM()',uniLM::corrupt_data()')")
  }
  invisible(validate_LM(data))
  check_1num(nboot, "num")

  formula<- momOverwrite("formula")
  beta<- momOverwrite("beta")
  data<-data$data

  # bootstrap logic:

  boots<-list(); length(boots) <- nboot

  if (nboot>0){
    for (i in 1:nboot){
      boots[[i]]<-calculate_mom(data=data, formula=formula, beta=beta, K=K,
                                algorithm=algorithm, stochastic=stochastic, ...)
    }

    return(list(
      betas=sapply(boots, function(x) x$b),
      b=rowMeans(betas),
      vcov = cov(t(betas)),
      se=sqrt(diag(vcov)),
      scores=rowMeans(sapply(boots,function(x)x$scores))
    ))

  } else{
    calculate_mom(data=data, formula=formula, beta=beta, K=K,
                  algorithm=algorithm, stochastic=stochastic, ...)
  }

}







