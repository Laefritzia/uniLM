
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
#' @param blocks supply block-partition of the data (eg fixed blocks or reshuffled each iteration)
#' @examples
#'  data<-uniLM::corrupt_data(n=100,scenario="a",O=0.1)
#'  X<-model.matrix(data$form,data=data$data)
#   y<-data$data[,all.vars(data$form)[1]])
#
#'  med_block(X=X,y=y,K=20,
#'  b=rep(0,length(data$beta)),
#'  b_prime=data$beta
#'  )
#' @returns named list with X and y of median block
#' @export
med_block <- function(X,y,K,b,b_prime, blocks=NULL){

  blocks<-if (is.null(blocks)) sample(factor(rep(1:K, length.out=nrow(X)))) else blocks

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
calculate_mom <- function(data, form, beta, K, algorithm=c("GD", "ADMM"),
                          stochastic=TRUE, stepsize=0.01, maxiter=100, seed=NULL){

  # checks and paramaters

  invisible(mapply(
    check_1num, list(K, stepsize, maxiter), c("int", "num", "int"))
  )

  alg<-match.arg(algorithm)
  X<-model.matrix(form,data=data)
  y<-data[,all.vars(form)[1]]

  if(!is.null(seed)) set.seed(seed)
  fix_blocks<-if(stochastic) NULL else sample(factor(rep(1:K, length.out=nrow(X))))

  # > Gradient Descent -----------------
  if(alg%in% "GD"){ # Gradient Descent

    b<-b_prime<-numeric(ncol(X))
    mom_obj<-mom_err<-numeric(maxiter)
    scores<-numeric(nrow(X))

    iter<-0
    while(TRUE) {
      iter<-iter+1

      blocks<- if(stochastic) sample(factor(rep(1:K, length.out=nrow(X)))) else fix_blocks
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
      mom_obj[iter]<- sum((X%*%b-y)^2)/nrow(X) - sum((X%*%b_prime-y)^2)/nrow(X)
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
    scores<-numeric(nrow(X))
    rho<-5 # same as in their paper?
    iter<-0
    while(TRUE) {
      iter<-iter+1
      blocks<- if(stochastic) sample(factor(rep(1:K, length.out=nrow(X)))) else fix_blocks
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
      mom_obj[iter]<- sum((X%*%b-y)^2)/nrow(X) - sum((X%*%b_prime-y)^2)/nrow(X)
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
#' @param form formula specifiyng the relationship of response and predictors in data
#' @param beta true coefficient for relationship (for simulation purposes)
#' @param K integer specifying number of blocks for MOM-algorithm
#' @param algorithm specifying how to get coefficients, GD (gradient-descent), ADMM(ascent-descent)
#' @param stochastic TRUE/FALSE decides, if blocks should be shuffled randomly in every iteration (stochastic==TRUE)
#'                   or based on a fixed starting partition (stochastic==FALSE)
#' @param ... optional arguments like stepsize, maxiter, seed, ...
#' @returns named list with final b (coefficients), iterative objectives and errors
#' @examples
#'   data<-uniLM::corrupt_data(n=100,scenario="a",O=0.1)
#'   res<-MOM(data, K=12, algorithm="GD")
#'   plot(0:100, res$mom_obj)
#'   plot(0:100, res$mom_err)
#'   abs(res$b-true_beta)
#' @export
MOM <- function(data, form=NULL, beta=NULL, K, algorithm=c("GD", "ADMM"),
                stochastic=TRUE, ...){
  UseMethod("MOM")
}


#' @rdname MOM
#' @export
MOM.default <- function(data, form, beta, K, algorithm=c("GD", "ADMM"),
                        stochastic=TRUE, ...){
  invisible(validate_LM(new_LM(list(data=data,form=form,beta=beta))))

  calculate_mom(data=data, form=form, beta=beta, K=K, algorithm=algorithm, stochastic=stochastic, ...)
}

#' @rdname MOM
#' @export
MOM.LM <- function(data, form=NULL, beta=NULL, K, algorithm=c("GD", "ADMM"),
                   stochastic=TRUE, ...){
  if(!inherits(data, "LM")){
    stop("Data must be of class LM (eg 'uniLM::LM()',uniLM::corrupt_data()')")
  }
  invisible(validate_LM(data))

  momOverwrite<-function(var){
    if(is.null(get(var,envir=parent.frame()))){
      return(data[[var]])
    } else{
      warning(paste0(var, " specified, overwriting value stored in LM-object (data)."))
    }
  }
  form<- momOverwrite("form")
  beta<- momOverwrite("beta")
  data<-data$data

  calculate_mom(data=data, form=form, beta=beta, K=K, algorithm=algorithm, stochastic=stochastic, ...)
  }







