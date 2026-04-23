
# > 1a ------------------
# internal constructor
new_LM <- function(res){
  obj <- structure(res,class="LM")
  obj
}

# > 1b ------------------
# for internal checks as well as converting formula and factor
validate_LM <- function(x){

  if(sum(names(x)%in%c("data","beta", "formula"))!=3){
    stop("Named list with at least data, beta, and formula must be supplied,\n
       create object of class 'LM', like with uniLM::corrupt_data() or uniLM::LM()")
  }

  if (!inherits(x$data, "data.frame")){
    stop("Data must be data.frame object.")
  }

  if(nrow(x$data)<2){
    stop("To few data points regression.")
  }


  ccols<-sapply(x$data,is.character)
  if (any(ccols)) {
    tryCatch({
      x$data[ccols] <- lapply(x$data[ccols],as.factor)
    }, error = function(e) {
      stop("Conversion of character to factor columns failed. Please supply\n
            numeric or factor columns only.")
    })
  }

  x$formula<-tryCatch(
    as.formula(x$formula),
    error=function(e){
      stop("Must supply regression formula in string/formula-format")
    }
  )

  if(any(!(all.vars(x$formula) %in% names(x$data)))){
    stop("Supplied formula does not match variables in data.")
  }

  if (!all(sapply(x$data,function(cols) is.numeric(cols) | is.factor(cols)))
      ) {
    stop("Only numeric or factor columns allowed.")
  }


  if(!is.null(x$beta) & !is.numeric(x$beta)){
    stop("Beta must be NULL or a numeric vector")
  }


frml_len <- ncol(model.matrix(x$formula, x$data[1,,drop=FALSE]))
 if(length(x$beta) != frml_len){
   stop("Regression formula doesnt match length of true beta.\nCheck if beta and formula match with your intercept logic and drop if not needed.")
 }

return(x)
}

# > 1c ------------------
#' print method for class LM
#
#' @param x object of class LM \(contains true beta and formula\)
#' @returns something super duper amazing
#' @export
print.LM <- function(x, ...){

  cat("\nCall:", deparse(x$formula),"\n\nData:\n")
  str(x$data)
  cat("\nbeta:\n")
  print(x$beta)

}

# > 1c ------------------
#' print method for class LMfit
#
#' @param x object of class LM \(contains true beta and formula\)
#' @returns something super duper amazing
#' @import stats
#' @export
print.LMfit <- function(x, ...){

  nboot<-attributes(x)$nboot
  if(nboot > 0)cat("\nResults for bootstrapped variance and point estimator with", nboot, "runs.\n")

  cat("\nCall: ", deparse(x$formula),"\n\nminmax MOM-result:\n\n")

  out <- as.data.frame(matrix(0, length(x$beta), 6))
  rownames(out) <- colnames(model.matrix(x$formula, x$data[1,,drop=FALSE]))
  colnames(out) <- c("Estimate", "Std.Error", "Wald-z", "p-value", "Lower", "Upper")

  out[,1] <- x$MOM$b
  c<-stats::qnorm(0.975)

  out[,-1] <- if(nboot > 0){
      cbind(x$MOM$se, x$MOM$b/x$MOM$se, 2*(1-stats::pnorm(abs(x$MOM$b/x$MOM$se))),
                             x$MOM$b-c*x$MOM$se, x$MOM$b+c*x$MOM$se)
  } else {
    rep(NA_real_,length(x$beta)*5)
  }

  print(out)
}

# > 2 ------------------
#' Robust linear model
#'
#' @param data data.frame with regression data
#' @param formula formula for linear model
#' @param beta for simulation purposes: true linear relationship between response and predictors
#' @param MOMalgorithm specifies algorithm to find minimax mom-estimator
#' @param stochastic logical, controls if blocks in MOM-algorithm should be reshuffled each iteration
#' @param nboot integer controlling how many times MOM-algorithm should be repeated on a resampled (with replacement) dataset.
#' @param ... additional arguments passed to MOM-algorithm, see ?MOM for details.
#' @returns named list of class LM
#' @import stats MASS
#' @export
LM <- function(data, formula, beta=NULL,
               MOMalgorithm=c("GD", "ADMM"),
               stochastic=TRUE,
               nboot=100,...){

  # bare LM
  bLM<-if(!inherits(data, "LM")){
    validate_LM(new_LM(list(data=data,formula=formula,beta=beta)))
  } else {
    validate_LM(data)
  }

  X <- model.matrix(bLM$formula, bLM$data[1,,drop=FALSE])

  if(is.null(data$beta)){
    data$beta<-rep(0, ncol(X))
  }

  if (nboot%in%0){
    warning("No variance for minmax-MOM calculated. Set nboot to an integer > 0 to get bootstrapped results.")
  }



  alg<-match.arg(MOMalgorithm)

  # I don't think we need to, but this is how we could force intercept:
  # mod<-terms(bLM$formula)
  # attributes(mod)$intercept <- 1
  # bLM$formula <- mod

  # choose best block-size for MOM
  K<-adaptK(bLM, algorithm=alg,stochastic=stochastic)

  # minimax mom method
  mod_MOM<-MOM(bLM, K=K, algorithm=alg, stochastic=stochastic, nboot=nboot,...)


  # estimators other than mom
  mod_lm<-stats::lm(bLM$formula, bLM$data)
  mod_rlm<-MASS::rlm(bLM$formula, bLM$data)

  res <- new_LM(
    list(
    MOM = mod_MOM,
    OLS = list(
      b=mod_lm$coeff,
      se=sqrt(diag(vcov(mod_lm)))
    ),
   M_est = list(
     b = mod_rlm$coefficients,
     se = sqrt(diag(vcov(mod_rlm)))
   ),
   # this is the input
   data = bLM$data,
   formula = bLM$formula,
   beta = bLM$beta
  )
  )

  # additional subclass LMfit
  class(res) <- c("LMfit", "LM")

  attr(res, "nboot") <- nboot

  return(res)

}

# > 2.2 adaptive K blocks -----------------

#' choose the best block-partition for MOM-algorithm in LM-framework
#'
#' @param data data.frame of class LM (in simplest form a list containing data, formula and beta)
#' @param formula formula for linear model
#' @param beta for simulation purposes: true linear relationship between response and predictors
#' @param K_grid integers specifying number of blocks for MOM-algorithm
#' @param algorithm specifying how to get coefficients, GD (gradient-descent), ADMM(ascent-descent)
#' @param stochastic TRUE/FALSE decides, if blocks should be shuffled randomly in every iteration (stochastic==TRUE)
#'                   or based on a fixed starting partition (stochastic==FALSE)
#' @param maxiter number of iterations the MOM-algorithm should perform
#' @returns an integer K specifiyng the ideal number of blocks to partition present data-structure.
#' @export
adaptK <- function(data, formula=NULL, beta=NULL,
                   K_grid = 2:(nrow(data$data)/2),
                   algorithm=c("GD", "ADMM"), maxiter=20,
                   stochastic=TRUE, ...){
  UseMethod("adaptK")
}


#' @rdname adaptK
#' @export
adaptK.default <- function(data, formula, beta, ...){

  # try to convert as LM object and pass to LM method:
  obj <- validate_LM(new_LM(list(data=data,formula=formula,beta=beta)))

  adaptK(obj)

}

#' @rdname adaptK
#' @export
adaptK.LM <- function(data, formula=NULL, beta=NULL,
                      K_grid = 2:(nrow(data$data)/2),
                      algorithm=c("GD", "ADMM"), maxiter=20,
                      stochastic=TRUE, ...){
  # Checks:
  if(!inherits(data, "LM")){
    stop("Data must be of class LM (eg 'uniLM::LM()',uniLM::corrupt_data()')")
  }
  data <- validate_LM(data)
  invisible(sapply(K_grid, check_1num, "int"))

  alg<-match.arg(algorithm)
  # grid search returning euclidean norm between estimate and beta

  grid<-sapply(K_grid, function(K){
    res<-MOM(data, K=K, algorithm=alg,
             maxiter=maxiter,stochastic=stochastic)

    sum((res$b-data$beta)^2)
  })

  # returning best K
  return(K_grid[which.min(grid)])

}



# > 2.3 ------------------
#' summary method for class LM
#'
#' @param data object of class LM (contains true beta and formula)
#' @returns something super duper amazing
#' @export
summary.LMfit <- function(data){

  # Checks:
  if(!inherits(data, "LMfit")){
    stop("Data must be of class LMfit (call 'uniLM::LM()'")
  }
  data<-validate_LM(data)

  # hier z und p rein?

}

# > 3 ------------------
#' residuals method for class LM
#'
#' @param data object of class LM \(contains true beta and formula\)
#' @param type string specifying if response, standardized, or studentized residuals
#' @returns vector of residuals
#' @examples
#'   data<-uniLM::corrupt_data(n=100,scenario="a",O=0.1)
#'   r_t <- residuals(data,type="student")
#'   head(r_t)
#' @export
residuals.LM <- function(data, type=c("response", "standard", "student")){

  # Checks:
  if(!inherits(data, "LM")){
    stop("Data must be of class LM (eg 'uniLM::LM()',uniLM::corrupt_data()')")
  }
  data<-validate_LM(data)
  t<-match.arg(type)

  # my objects
  n<-nrow(data$data)
  p<-length(data$beta)-1
  y<-data$data[,all.vars(data$formula)[1]]
  X<-model.matrix(data$formula, data=data$data)
  # algebra
  H<-X%*%solve(t(X)%*%X)%*%t(X)   # projectionmatrix on span(X)
  M<- diag(n)-H                   # orthogonal complement for H
  ols<-solve(t(X)%*%X)%*%t(X)%*%y # ordinary-least-squares-estimator
  yhat<-X%*%ols # H%*%y

  # response residuals
  ehat<-y-yhat
  # standardized residuals
  shat<-sum(ehat^2)/(n-(p+1))
  r_std<-ehat/(sqrt(shat*(1-diag(H))))
  # studentized residuals (umwandlung: https://online.stat.psu.edu/stat462/node/247/)
  r_t<- r_std*sqrt((n-p-2)/(n-p-1-r_std^2))

  return(switch(t,
         "response"=ehat,
         "standard"=r_std,
         "student" =r_t))

}

# > 4 ------------------
#' diagnostic plots for class LM
#'
#' @param data object of class LM \(contains true beta and formula\)
#' @param plot TRUE draws plot, FALSE stores diagnostic values
#' @param which can specify plots to be drawn
#' @returns plot output or named list with values
#' @examples
#'   data<-uniLM::corrupt_data(n=100, scenario="a",O=0.1)
#'   plot(data,which=3)
#'
#'   plot_data<-plot(data,plot=FALSE)
#'   head(plot_data$cook_d)
#'
#' @import ggplot2 stats
#' @export
plot.LM <- function(data, plot=TRUE, which=1:4, ...){

  # Checks
  if(!inherits(data, "LM")){
    stop("Data must be of class LM (eg 'uniLM::LM()',uniLM::corrupt_data()')")
  }
  data<-validate_LM(data)
  if(any(!(which %in% 1:4))) stop("If 'which'-argument, must specify 1-4 for type of plot.")

  # Data
  ehat<-residuals(data, type="response")
  y<-data$data[,all.vars(data$formula)[1]]
  X<-model.matrix(data$formula, data=data$data)
  n<-nrow(data$data)
  p<-length(data$beta)-1
  H<-X%*%solve(t(X)%*%X)%*%t(X)
  yhat<-y-ehat
  r_std<-residuals(data, type="standard")
  r_t  <-residuals(data, type="student")
  # Cooks Distance (https://de.wikipedia.org/wiki/Cook-Abstand)
  cook_d<-r_t^2/(p+1)*(diag(H)/(1-diag(H)))


  if(plot){

    # Plot 1: fitted vs standardized studentized res
    g1<- data.frame(x=sqrt(abs(r_t)), y=yhat) |>
      ggplot2::ggplot(ggplot2::aes(x, y))+
      ggplot2::geom_point(alpha=0.5) +
      ggplot2::geom_line(data=data.frame(
        x=sqrt(abs(r_t)),
        y=stats::lowess(sqrt(abs(r_t)), yhat)$y),
        ggplot2::aes(x,y),color="Darkblue") +
      ggplot2::labs(x="sqrt(abs(studentized residuals))",
           title="fitted~sqrt(abs(r_t))")+
      ggplot2::theme_minimal()


    # Plot 2: studentized res vs. leverage
    g2 <- data.frame(x=diag(H), y=r_t) |>
      ggplot2::ggplot(ggplot2::aes(x,y))+
      ggplot2::geom_point(alpha=0.5) +
      ggplot2::labs(title="Studentized residuals~leverage",
           x="diag(H)", y="r_t")+
      ggplot2::theme_minimal()


    # Plot 3: Cook plot (red line is treshhold for outlier)
    g3 <- data.frame(x=1:n, y=cook_d) |>
      ggplot2::ggplot(ggplot2::aes(x,y))+
      ggplot2::geom_col() +
      ggplot2::geom_hline(yintercept=4/n, color="red") +
      ggplot2::geom_text(
        ggplot2::aes(label=ifelse(y>4/n,as.character(x),"")),
        vjust=-0.8,size=3.5,color = "black"
      ) +
      ggplot2::ylim(c(0, max((4/n)*1.02,max(cook_d)*1.02)))+
      ggplot2::labs(x="observations (row index)", y="Cooks Distance",
           title="Cooks Distance",
           subtitle=paste0("outlier idx: ",paste(as.character(1:n)[cook_d>4/n],collapse=","))
      )+
      ggplot2::theme_minimal()

    # Plot 4: QQ-Plot: empirical vs theoretical distribution
    yq <- quantile(ehat,c(0.25, 0.75))
    xq <- qnorm(c(0.25, 0.75))
    a <- diff(yq)/diff(xq)
    b <- yq[1]-a*xq[1] # y= a*x + b

    # vielleicht kann man das noch besser zoomen
    g4 <- data.frame(x=qnorm((1:n)/n - 0.01), y=sort(ehat)) |>
      ggplot2::ggplot(ggplot2::aes(x,y)) +
      ggplot2::geom_point(alpha=0.5)+
      ggplot2::geom_abline(ggplot2::aes(intercept=b,slope=a),color="Darkblue")+#color="robust"))+
      #ggplot2::geom_abline(ggplot2::aes(intercept=mean(ehat),slope=sd(ehat),color="naive"))+
      ggplot2::labs(x="theoretical quantile", y="residuals", #colour="",
           title="QQ Plot",
           subtitle="empirical~theoretical normal distribution"
      ) +
      ggplot2::theme(legend.position="bottom")+
      ggplot2::theme_minimal()

    plots<-list(g1,g2,g3,g4)

    for(i in which) print(plots[[i]])

  }

  return(
    invisible(cbind(data.frame(
   idx=1:n,
   y=y, yhat=yhat, ehat=ehat,
   leverage=diag(H),
   ehat=ehat, r_std=r_std, r_t=r_t,
   cook_d=cook_d), X))
  )

}

# > 5 ------------------
#' @export
confint.LMfit <- function(object, parm = "", level = 0.95,
                           test=c("", ""),...){


      test <- match.arg(test)
}

# > 6 -------------------------------

#' @export
as.data.frame.LM <- function(x, ...) {
  return(as.data.frame(x$data, ...))
}



