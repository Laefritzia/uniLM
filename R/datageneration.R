# -Scenario A: Constructed completely corrupted datapoints that directly contradict the true model (negative slope,..)
# -Scenario B: heavy tail MODEL und random outlier
# -Scenario C: perfektes model und random outlier
# -Scenario D: perfektes model nur mit random noise
#
# Scenario B:
# Goerg, Georg M., The Lambert Way to Gaussianize Heavy-Tailed Data with
# the Inverse of Tukey’s h Transformation as a Special Case,
# The Scientific World Journal, 2015, 909231, 16 pages, 2015.
# https://doi.org/10.1155/2015/909231
#
# -Trying cauchy errors first (rnorm/rnorm)

#' generates \(corrupt\) data of class LM for different scenarios
#'
#' @param n integer specifiyng number of observations
#' @param scenario as specified in seminar-paper
#' @param beta true linear relationship between response and predictors for CLEAN data (Inliers)
#' @param O in 0,1: percentage of data to be corrupted (O..Outliers)
#' @param OO worst range for Outliers
#' @param eps_var variance of random errors
#' @returns named list of class LM
#' @examples
#'   data<-corrupt_data(n=100, scenario="a",O=0.1)
#'   head(data$data)
#'   data$beta
#'   data$form
#'
#' @importFrom stats rnorm
#' @export
corrupt_data <- function(n, scenario=c("a", "b", "c", "d"),
                         beta=c(1,3),
                         O=0.2, OO=100, eps_var=0.3
){

  # a)
  # clean data (Inliers) follows true beta + gaussian noise
  # corrupt data (Outliers: cauchy and random): shifted intercept, opposite slope
  # b)
  # clean data (Inliers) follows true beta + cauchy noise to be heavy tailed
  # corrupt data (Outliers: cauchy and random) follow no particular distribution( and relationship with y?)
  #c)
  # clean data (Inliers) follows true beta + random noise
  # corrupt data (Outliers: cauchy and random) follow no particular distribution( and relationship with y?)
  # d)
  # all data follows true beta + random noise



  # Checks
  invisible(
    mapply(check_1num, list(n,O,OO,eps_var),c("num", "prob", rep("pos_num",2)))
  )

  scenario<-match.arg(scenario)

  eps<-rnorm(n)*eps_var
  eps_cauchy<- (rnorm(n)/rnorm(n)) *eps_var

  # drawing (and shuffling) Inliers and Outliers based on O-probabilities
  pos<-factor(
    sample(1:3, n, prob = c(1-O, O/2, O/2), replace=TRUE),
    levels=1:3,
    labels=c("I", "O1", "O2")
  )

  # generating corrupted vector
  x2<-numeric(n)
  x2_good <- rnorm(n)
  x2[pos%in%"I"]<-rnorm(sum(pos %in% "I")) #normal
  x2[pos%in%"O1"]<-rnorm(sum(pos %in% "O1"))/rnorm(sum(pos %in% "O1")) #cauchy
  x2[pos%in%"O2"]<-sample(seq(-OO,OO,by=0.001), sum(pos %in% "O2"), replace=TRUE) #random

  # building the model
  formula<-as.formula("y~x2") # formula
  clean<-as.integer(pos%in%"I") # index clean data
  X <- model.matrix(formula, data=data.frame(x2=x2,y=numeric(n))) # Intercept+covariates
  X_good <- model.matrix(formula, data=data.frame(x2=x2_good,y=numeric(n)))

  # defining the scenario
  y <- switch(scenario,
              "a" = clean*(X%*%beta+eps) + (1-clean)*(OO+beta[1]+X[,-1,drop=FALSE]%*%-beta[2]+eps),
              "b" = X%*%beta+eps_cauchy,
              "c" = X%*%beta+eps,
              "d" = X_good%*%beta+eps
  )

  x2 <- if(scenario %in% "d") x2_good else x2

  # returning LM object
  res<-new_LM(
    list(
      data=data.frame(x2=x2,eps=eps,y=y, clean=clean),
      beta=beta,
      formula=formula
    ))

  return(validate_LM(res))
}





# small ball :)
