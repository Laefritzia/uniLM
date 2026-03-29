
# > 1a ------------------
new_LM <- function(res){# internal constructor
  obj <- structure(res, class="LM")
  #attributes(obj, "strat_var") <- strat_var # kann hier noch Sachen reingebn
  obj
}

# > 1b ------------------
validate_LM <- function(x){ # for internal checks
  if (!inherits(x$data, "data.frame")){
    stop("Data must be data.frame object.")
  }
  if(any(!sapply(x$data,is.numeric))) {
    stop("Data must only consist of numeric columns. Use factor() for ranked/binary data.")
  }
  if(nrow(x$data)<2){
    stop("To few data points regression.")
  }

  if(sum(names(x)%in%c("data","beta", "form"))!=3){
  stop("Named list with data, true beta, and formula must be supplied,\n
       create object of class 'LM', like with uniLM::corrupt_data() or uniLM::LM()")
  }

  if(!is.numeric(x$beta)){
    stop("Beta must be a numeric vector")
  }

 if(length(x$beta) != length(all.vars(x$form))){
   stop("Regression formula doesnt match length of true beta.\nBeta must include intercept???")
 }
x
}

# > 2 ------------------
#' Robust linear model
#'
#' @param form formula for linear model
#' @param data data.frame with regression data
#' @param beta for simulation purposes: true linear relationship between response and predictors
#' @returns named list of class LM
#' @export
LM <- function(form, data, beta=NULL){
  # To-DO:
  # output und form des outputs reingeben
  # sollten mit intercept dealen (attributes(terms(formula)))

  form<-tryCatch(
    as.formula(form),
    error=function(e){
      stop("Must supply regression formula in string/formula-format")
    }
    )

  if(is.null(beta)){
    beta<-rep(0, length(all.vars(form))) # deal with intercept
  }

  res<-new_LM(
    list(data=data,beta=beta,form=form)
  )

  # wollen returnen: coefficient (MOM), resiudals, fitted,
  # varianz: welche?
  # test statistic und confidence interval: welche?

  # davor hier dann noch output definieren:
  return(validate_LM(res))

}

# > 3 ------------------
#' residuals method for class LM
#'
#' @param data object of class LM \(contains true beta and formula\)
#' @param type string specifying if response, standardized, or studentized residuals
#' @returns vector of residuals
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
  y<-data$data[,all.vars(data$form)[1]]
  X<-model.matrix(data$form, data=data$data)
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
         "student" =r_std))

}

# > 4 ------------------
#' diagnostic plots for class LM
#'
#' @param data object of class LM \(contains true beta and formula\)
#' @param plot TRUE draws plot, FALSE stores diagnostic values
#' @returns plot output or named list with values
#' @import ggplot2 tibble stats
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
  y<-data$data[,all.vars(data$form)[1]]
  X<-model.matrix(data$form, data=data$data)
  n<-nrow(data$data)
  p<-length(data$beta)-1
  H<-X%*%solve(t(X)%*%X)%*%t(X)
  yhat<-y-ehat
  r_std<-residuals(data, type="standard")
  r_t  <-residuals(data, type="student")
  # Cooks Distance (https://de.wikipedia.org/wiki/Cook-Abstand)
  cook_d<-r_t^2/(p+1)*(diag(H)/(1-diag(H)))


  if(plot){

    p_out<-which %in% c(1:4)

    # Plot 1: fitted vs standardized studentized res
    g1<- tibble(x=sqrt(abs(r_t)), y=yhat) |>
      ggplot(aes(x, y))+
      geom_point(alpha=0.5) +
      geom_line(data=tibble(
        x=sqrt(abs(r_t)),
        y=stats::lowess(sqrt(abs(r_t)), yhat)$y),
        aes(x,y),color="Darkblue") +
      labs(x="sqrt(abs(studentized residuals))",
           title="fitted~sqrt(abs(r_t))")


    # Plot 2: studentized res vs. leverage
    g2 <- tibble(x=diag(H), y=r_t) |>
      ggplot(aes(x,y))+
      geom_point(alpha=0.5) +
      labs(title="Studentized residuals~leverage",
           x="diag(H)", y="r_t")


    # Plot 3: Cook plot (red line is treshhold for outlier)
    g3 <- tibble(x=1:n, y=cook_d) |>
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
    g4 <- tibble(x=qnorm((1:n)/n - 0.01), y=sort(ehat)) |>
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
# basierend auf? irgendwas robustem
confint.LM <- function(object, data, parm = "trt", level = 0.95,
                           test=c("Score", "Wald"),
                           tol=10^-9, maxiter=20, ...){

      if (!inherits(object, "my_cox")){
        stop("Object must be class 'my_cox'.")
      }

      test <- match.arg(test)
}






