
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
  med_ind<-which.min(abs(stats::median(means_loss, na.rm=TRUE)-means_loss))

  if(is.na(med_ind)||length(med_ind)==0) med_ind <- 1

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

  mf<-model.frame(formula,data)
  y<-mf[,1]
  X<-model.matrix(formula, data)
  n<-nrow(X)

  if(!is.null(seed)) set.seed(seed)
  fix_blocks<-if(stochastic) NULL else sample(factor(rep(1:K, length.out=n)))


  b<-b_prime<-u<-u_prime<-z<-z_prime<-numeric(ncol(X))
  mom_obj<-mom_err<-numeric(maxiter)
  scores<-numeric(n)

  # > Gradient Descent -----------------
  if(alg%in% "GD"){ # Gradient Descent

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
  }

  # > ADMM -----------------
  if(alg%in% "ADMM"){

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

}

  return(
      structure(
        list(
          b=b,
          mom_obj=mom_obj,
          mom_err=mom_err,
          scores=scores/maxiter,
          maxiter=maxiter,
          K=K,
          n=n,
          se=rep(NA_real_, ncol(X))
        ),
        class="MOM"
      )
    )

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
#' @param nboot integer controlling how many times MOM-algorithm should be repeated on a resampled (with replacement) dataset.
#' @param parallel TRUE/FALSE decides, if the bootstrapping (nboot>0) should be parallelized. Only use this if you know how to plan sequential vs. multisession
#' @param ... optional arguments like stepsize, maxiter, seed, ...
#' @returns named list with final b (coefficients), iterative objectives and errors
#' @import utils future.apply
#' @export
MOM <- function(data, formula=NULL, beta=NULL, K, algorithm=c("GD", "ADMM"),
                stochastic=TRUE, nboot=0, parallel=FALSE,...){
  UseMethod("MOM")
}


#' @rdname MOM
#' @export
MOM.default <- function(data, formula, beta, ...){

  # try to convert as LM object and pass to LM method:
  obj<-validate_LM(new_LM(list(data=data,formula=formula,beta=beta)))

  MOM(obj, ...)

}

#' @rdname MOM
#' @export
MOM.LM <- function(data, formula=NULL, beta=NULL, K, algorithm=c("GD", "ADMM"),
                   stochastic=TRUE, nboot=0, parallel=FALSE,...){

  # checks and structure:

  if(!inherits(data, "LM")){
    stop("Data must be of class LM (eg 'uniLM::LM()',uniLM::corrupt_data()')")
  }
  data<-validate_LM(data)
  invisible(check_1num(nboot, "num"))

  formula<- if(is.null(formula)) data$formula else warning("External formula object, other than the one stored in LM-object, will be used")
  beta<- if(is.null(beta)){
    data$beta
    } else {
      message("Beta in LM-object will now contain the minmax MOM-estimate, starting the algorithm at: ", beta)
    }

  data<-data$data


  RES<-calculate_mom(data=data, formula=formula, beta=beta, K=K,
                     algorithm=algorithm, stochastic=stochastic, ...)

  # 'bootstrap-like' resampling, but but just reshuffle blocks on same dataset (effectively an empirical variance)

  boots<-list(); length(boots) <- nboot

  if (nboot>0){

    # only give the pretty message, when not parallelizing
    if (!parallel){
      message("Bootstrapping progress: \n")

    for (i in 1:nboot){
      data_b<- data[sample(1:nrow(data),replace=TRUE) ,]
      boots[[i]]<-calculate_mom(data=data_b, formula=formula, beta=beta, K=K,
                                algorithm=algorithm, stochastic=stochastic, ...)

      #if(i%%5 == 0){
      cat("\r",round(i/nboot*100), "%", sep="") #\r always restarts on the same line
      utils::flush.console()
      #}
    } } else {
      boots<-future.apply::future_lapply(1:nboot, function(x){
        data_b<- data[sample(1:nrow(data),replace=TRUE) ,]
        calculate_mom(data=data_b, formula=formula, beta=beta, K=K,
                      algorithm=algorithm, stochastic=stochastic, ...)
      }, future.seed=TRUE)

    }

    betas<-sapply(boots, function(x) x$b)
    beta_ok <- colSums(is.finite(betas))==nrow(betas) # one row in beta is one coeff
    vcov_b <- cov(t(betas[,beta_ok,drop=FALSE]))

    return(
      structure(list(
      b = RES$b,
      mom_obj=RES$mom_obj,
      mom_err=RES$mom_err,
      scores=RES$scores/RES$maxiter,
      maxiter=RES$maxiter,
      K=RES$K,
      n=RES$n,
      betas=betas,
      boot_b=rowMeans(betas),
      vcov = vcov_b,
      se = sqrt(diag(vcov_b)),
      nboot=nboot,
      nboot_conv = sum(beta_ok)
    ),
    class="MOM")
    )

  } else{

    RES

  }

}



# > plot MOM --------------------------------------

#' plot results from MOM-algorithm, a.o. useful for outlier detection
#'
#' @param mom object containing results from MOM-algorithm
#' @param plot TRUE draws plot, FALSE stores results (including outlier detection)
#' @param which can specify plots to be drawn
#' @param block_p upper bound for scores, everything below is classified as outlier
#' @param ... additional arguments passed to plot (currently not used)
#' @returns plot output or named list with value
#' @importFrom rlang .data
#' @import ggplot2
#' @export
plot.MOM <- function(mom, plot=TRUE, which=1:4, block_p=0.001, ...){

  algorithm_data <- data.frame(
    iterations=1:mom$maxiter,
    mom_err=mom$mom_err,
    mom_obj=mom$mom_obj
  )

  outlier_data <-data.frame(
    x=1:mom$n, scores=mom$scores
    )

  #block_p <- 0.001
  #(1/mom$K) - sqrt((1-1/mom$K)*(1/mom$K)/mom$n) #minus expected SE (variance of a frequncy: binomial)


  outlier_data$outlier <- factor(as.integer(outlier_data$scores<=block_p),
                              levels=c(0,1), labels=c("Inlier", "Outlier")
  )

  ylimit<-max(block_p, outlier_data$scores)*1.2


  if(plot){

suppressMessages({suppressWarnings({
    g1<- algorithm_data |>
      ggplot2::ggplot(ggplot2::aes(.data$iterations, .data$mom_err))+
      ggplot2::geom_line(color="Skyblue",lwd=1.5)+
      ggplot2::labs(title="Error between estimate and true")+
      ggplot2::theme_minimal()

    g2<- algorithm_data |>
      ggplot2::ggplot(ggplot2::aes(.data$iterations, .data$mom_obj))+
      ggplot2::geom_line(color="Skyblue",lwd=1.5)+
      ggplot2::labs(title="MOM-objective (distance between the two candidates)")+
      ggplot2::theme_minimal()

    g3 <- outlier_data |>
      ggplot2::ggplot(ggplot2::aes(.data$x, .data$scores, fill=.data$scores<=block_p)) +
      ggplot2::geom_col() +
      ggplot2::geom_point(ggplot2::aes(color=.data$scores<=block_p))+
      ggplot2::geom_hline(yintercept=block_p, color="red")+
      ggplot2::scale_fill_manual(values=c(
        "TRUE"="#be0032",
        "FALSE"="#007bb8"
      ))+
      ggplot2::scale_color_manual(values=c(
        "TRUE"="#be0032",
        "FALSE"="#007bb8"
      ))+
      ggplot2::geom_text(
        ggplot2::aes(label=ifelse(.data$scores<=block_p, as.character(.data$x), "")),
        vjust=-0.8,size=3.5,color="black"
      ) +
      ggplot2::ylim(c(0, ylimit))+
      ggplot2::labs(title="Outlier Detection via Median-Block-Frequency",
           subtitle=paste0("outlier idx: ",paste(as.character(1:mom$n)[outlier_data$scores<=block_p],collapse=","))
           )+
       ggplot2::theme_minimal()+
      ggplot2::theme(legend.position = "bottom")

    g4 <- outlier_data |>
      ggplot2::ggplot(ggplot2::aes(.data$outlier,.data$scores, fill=.data$outlier)) +
      ggplot2::geom_boxplot() +
      #ggplot2::geom_hline(yintercept=block_p, color = "red") +
      ggplot2::labs(title = "Scores grouped by Outlier Status",
           subtitle="Scores by Median-Block-Frequency",x = "Status", y = "scores",
           fill="")+
      ggplot2::theme_minimal()

    plots<-list(g1,g2,g3,g4)

})})

    for(i in which) print(plots[[i]])

  }

  return(
    invisible(list(algorithm_data=algorithm_data,
                   outlier_data = outlier_data))
  )


}



