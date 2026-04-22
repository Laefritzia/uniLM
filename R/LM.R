
# > 1a ------------------
# internal constructor
new_LM <- function(res){
  obj <- structure(res, class="LM")
  #attributes(obj, "strat_var") <- strat_var # kann hier noch Sachen reingebn
  obj
}

# > 1b ------------------
# for internal checks
validate_LM <- function(x){
  if (!inherits(x$data, "data.frame")){
    stop("Data must be data.frame object.")
  }
  if(any(!sapply(x$data,is.numeric))) {
    stop("Data must only consist of numeric columns. Use factor() for ranked/binary data.")
  }

  # x$formula<-tryCatch(
  #   as.formula(formula),
  #   error=function(e){
  #     stop("Must supply regression formula in string/formula-format")
  #   }
  # )

  # nfCols_bf<-sapply(x$data, function(d) is.numeric(d) | is.factor(d))
  # x$data[,!nfCols] <- tryCatch(
  #   { cCols<-names(x$data[,!nfCols,drop=FALSE])
  #    sapply(cCols, function(c){
  #      factor(x$data[[c]])
  #    })
  #   },
  #   error=function(e){
  #     stop("Conversion of character- to factor-columns failed. Please supply\n
  #          numeric or factor cols")
  #   }
  # )
  # nfCols_af<-sapply(x$data, function(d) is.numeric(d) | is.factor(d))
#
#   if(any(!nfCols_ad)){
#     stop("Data contains invalid columns (character, numeric or factor allowed).")
#   }

  if(nrow(x$data)<2){
    stop("To few data points regression.")
  }

  if(sum(names(x)%in%c("data","beta", "formula"))!=3){
  stop("Named list with data, true beta, and formula must be supplied,\n
       create object of class 'LM', like with uniLM::corrupt_data() or uniLM::LM()")
  }

  if(!is.numeric(x$beta)){
    stop("Beta must be a numeric vector")
  }
#
#   x$formula<-tryCatch(
#     as.formula(formula),
#     error=function(e){
#       stop("Must supply regression formula in string/formula-format")
#     }
#   )


 if(length(x$beta) != length(all.vars(x$formula))){
   stop("Regression formula doesnt match length of true beta.\nBeta must include intercept???")
 }
x
}

# > 2 ------------------
#' Robust linear model
#'
#' @param formula formula for linear model
#' @param data data.frame with regression data
#' @param MOMalgorithm specifies algorithm to find minimax mom-estimator
#' @param beta for simulation purposes: true linear relationship between response and predictors
#' @returns named list of class LM
#' @import stats MASS
#' @export
LM <- function(data, formula, MOMalgorithm=c("GD", "ADMM"),stochastic=stochastic,
               nboot=0,beta=NULL){

  if(is.null(bLM$beta)){
    bLM$beta<-rep(0, ncol(X))
  }

  if (nboot%in%0){
    warning("Set nboot to any integer > 0 to get a (bootstrapped) variance for minimax MOM-estimator.")
  }

  # bare LM
  bLM<-validate_LM(new_LM(list(data=data,formula=formula,beta=beta)))

  # I don't think we need to, but this is how we could force intercept:
  # mod<-terms(bLM$formula)
  # attributes(mod)$intercept <- 1
  # bLM$formula <- mod

  # choose best block-size for MOM
  K<-adaptK(bLM, algorithm=match.arg(MOMalgorithm),stochastic=stochastic)

  # minimax mom method
  mod_MOM<-MOM(bLM, K=K, algorithm=MOMalgorithm, stochastic=stochastic, nboot=nboot)


  # estimators other than mom
  mod_lm<-stats::lm(bLM$formula, bLM$data)
  mod_rlm<-MASS::rlm(bLM$formula, bLM$data)

  res <- list(
    MOM = mod_lm,
    OLS = list(
      b=mod_lm$coeff,
      se=sqrt(vcov(mod_lm))
    ),
   M_est = list(
     b = mod_rlm$coefficients,
     se = sqrt(vcov(mod_rlm))
   )
  )

  # hier reingeben:
  # die drei estimators

  # robuste varianz: darf ich diese random robuste varianz vom MOM nehmen?

  # wollen returnen: coefficient (MOM), resiudals, fitted,
  # varianz: welche?
  # test statistic und confidence interval: welche?

  # result spalte reingeben
  # schoener output?
  return(validate_LM(new_LM(res)))

}

# > 2.2 adaptive K blocks -----------------

#' choose the best block-partition for MOM-algorithm in LM-framework
#'
#' @param data data.frame of class LM (in simplest form a list containing data, formula and beta)
#' @param K_grid integers specifying number of blocks for MOM-algorithm
#' @param algorithm specifying how to get coefficients, GD (gradient-descent), ADMM(ascent-descent)
#' @param stochastic TRUE/FALSE decides, if blocks should be shuffled randomly in every iteration (stochastic==TRUE)
#'                   or based on a fixed starting partition (stochastic==FALSE)
#' @param maxiter number of iterations the MOM-algorithm should perform
#' @returns an integer K specifiyng the ideal number of blocks to partition present data-structure.
#' @export
adaptK <- function(data,...){
  UseMethod("adaptK")
}

#' @rdname adaptK
#' @export
adaptK.LM <- function(data, K_grid = 2:(nrow(data$data)/2),
                      algorithm=c("GD", "ADMM"), maxiter=20,
                      stochastic=TRUE){
  # Checks:
  if(!inherits(data, "LM")){
    stop("Data must be of class LM (eg 'uniLM::LM()',uniLM::corrupt_data()')")
  }
  invisible(validate_LM(data))
  sapply(K_grid, check_1num, "int")

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

# > 2.2 ------------------
#' print method for class LM
#'
#' @param x object of class LM \(contains true beta and formula\)
#' @returns something super duper amazing
#' @export

# print.LM <- function(x, ...){
#
# }




# > 2.3 ------------------
#' summary method for class LM
#'
#' @param data object of class LM \(contains true beta and formula\)
#' @returns something super duper amazing
#' @export

summary.LM <- function(data){

  # Checks:
  if(!inherits(data, "LM")){
    stop("Data must be of class LM (eg 'uniLM::LM()',uniLM::corrupt_data()')")
  }
  invisible(validate_LM(data))

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
  invisible(validate_LM(data))
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

  args<-list(...) # keine ahnung was ich jetzt damit mach

  # Checks
  if(!inherits(data, "LM")){
    stop("Data must be of class LM (eg 'uniLM::LM()',uniLM::corrupt_data()')")
  }
  invisible(validate_LM(data))
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
      ggplot(aes(x, y))+
      geom_point(alpha=0.5) +
      geom_line(data=data.frame(
        x=sqrt(abs(r_t)),
        y=stats::lowess(sqrt(abs(r_t)), yhat)$y),
        aes(x,y),color="Darkblue") +
      labs(x="sqrt(abs(studentized residuals))",
           title="fitted~sqrt(abs(r_t))")


    # Plot 2: studentized res vs. leverage
    g2 <- data.frame(x=diag(H), y=r_t) |>
      ggplot(aes(x,y))+
      geom_point(alpha=0.5) +
      labs(title="Studentized residuals~leverage",
           x="diag(H)", y="r_t")


    # Plot 3: Cook plot (red line is treshhold for outlier)
    g3 <- data.frame(x=1:n, y=cook_d) |>
      ggplot(aes(x,y))+
      geom_col() +
      geom_hline(yintercept=4/n, color="red") +
      geom_text(
        aes(label=ifelse(y>4/n,as.character(x),"")),
        vjust=-0.8,size=3.5,color = "black"
      ) +
      ylim(c(0, max((4/n)*1.02,max(cook_d)*1.02)))+
      labs(x="observations (row index)", y="Cooks Distance",
           title="Cooks Distance",
           subtitle=paste0("outlier idx: ",paste(as.character(1:n)[cook_d>4/n],collapse=","))
      )

    # Plot 4: QQ-Plot: empirical vs theoretical distribution
    yq <- quantile(ehat,c(0.25, 0.75))
    xq <- qnorm(c(0.25, 0.75))
    a <- diff(yq)/diff(xq)
    b <- yq[1]-a*xq[1] # y= a*x + b

    # vielleicht kann man das noch besser zoomen
    g4 <- data.frame(x=qnorm((1:n)/n - 0.01), y=sort(ehat)) |>
      ggplot(aes(x,y)) +
      geom_point(alpha=0.5)+
      geom_abline(aes(intercept=b,slope=a),color="Darkblue")+#color="robust"))+
      #geom_abline(aes(intercept=mean(ehat),slope=sd(ehat),color="naive"))+
      labs(x="theoretical quantile", y="residuals", #colour="",
           title="QQ Plot",
           subtitle="empirical~theoretical normal distribution"
      ) +
      theme(legend.position="bottom")

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
# basierend auf zB sandwich vhcov? irgendwas robustem
confint.LM <- function(object, data, parm = "trt", level = 0.95,
                           test=c("Score", "Wald"),
                           tol=10^-9, maxiter=20, ...){

      if (!inherits(object, "my_cox")){
        stop("Object must be class 'my_cox'.")
      }

      test <- match.arg(test)
}

# > 6 -------------------------------
as.data.frame.LM <- function(x, ...) {
  # sonst ist data.frame Methode verwirrt von $beta und $formula
  return(as.data.frame(x$data, ...))
}



