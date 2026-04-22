# Scenario A: Constructed completely corrupted datapoints that directly contradict the true model (negative slope,..)
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
#' @param beta true linear relationship between response and predictors for CLEAN data (Inliers)
#' @param O in 0,1: percentage of data to be corrupted \(O..Outliers\)
#' @param OO worst range for Outliers
#' @param scenario as specified in seminar-paper
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
                         beta=list(a=c(1,3), b=NULL, c=NULL, d=NULL),
                         O=0.2, OO=100, eps_var=0.3
){

  # To-Do:
  # .) more scenarios: good data, heavy tails, multivariable,...
  # .) deal with intercept

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
    labels=c("I", "O1", "O2")
  )

  # generating corrupted vector
  x2<-numeric(n)
  x2[pos%in%"I"]<-rnorm(sum(pos %in% "I")) #normal
  x2[pos%in%"O1"]<-rnorm(sum(pos %in% "O1"))/rnorm(sum(pos %in% "O1")) #cauchy
  x2[pos%in%"O2"]<-sample(seq(-OO,OO,by=0.001), sum(pos %in% "O2"), replace=TRUE) #random

  if(scenario %in% "a"){
    # clean data (Inliers) follows true beta + gaussian noise
    # corrupt data (Outliers: cauchy and random): shifted intercept, opposite slope

    # building the model
    betaA<-beta$a
    formula<-as.formula("y~x2") # formula
    clean<-as.integer(pos%in%"I") # index clean data
    X <- model.matrix(formula, data=data.frame(x2=x2,y=numeric(n))) # Intercept+covariates

    y<- clean*(X%*%betaA+eps) + # Clean data
      (1-clean)*(OO+betaA[1]+X[,-1,drop=FALSE]%*%-betaA[-1]+eps) # corrupt data

    res<-new_LM(list(
      data=data.frame(x2=x2,eps=eps,y=y, clean=clean),
      beta=betaA,
      formula=formula))
    0
    return(validate_LM(res))
  }

  if(scenario %in% "b"){
    # clean data (Inliers) follows true beta + cauchy noise to be heavy tailed
    # corrupt data (Outliers: cauchy and random) follow no particular distribution( and relationship with y?)

    # building the model
    betaB<-beta$b
    formula<-as.formula("y~x2") # formula
    clean<-as.integer(pos%in%"I") # index clean data
    X <- model.matrix(formula, data=data.frame(x2=x2,y=numeric(n))) # Intercept+covariates

    y<- clean*(X%*%betaB+eps_cauchy) + # Clean data
      (1-clean)*X # corrupt data

res<-new_LM(list(
  data=data.frame(x2=x2,eps=eps,y=y, clean=clean),
  beta=betaB,
  formula=formula))

return(validate_LM(res))
}


  if(scenario %in% "c"){

    # clean data (Inliers) follows true beta + random noise
    # corrupt data (Outliers: cauchy and random) follow no particular distribution( and relationship with y?)

    # building the model
    betaB<-beta$C
    formula<-as.formula("y~x2") # formula
    clean<-as.integer(pos%in%"I") # index clean data
    X <- model.matrix(formula, data=data.frame(x2=x2,y=numeric(n))) # Intercept+covariates

    y<- clean*(X%*%betaC+eps) + # Clean data
      (1-clean)*X) # corrupt data

res<-new_LM(list(
  data=data.frame(x2=x2,eps=eps,y=y, clean=clean),
  beta=betaC,
  formula=formula))

return(validate_LM(res))
  }

  if(scenario %in% "d"){
    # all data follows true beta + random noise

    # building the model
    betaD<-beta$d
    formula<-as.formula("y~x2") # formula
    #clean<-as.integer(pos%in%"I") # index clean data
    X <- model.matrix(formula, data=data.frame(x2=x2,y=numeric(n))) # Intercept+covariates

    y<- X%*%betaD+eps

    res<-new_LM(list(
      data=data.frame(x2=x2,eps=eps,y=y, clean=clean),
      beta=betaD,
      formula=formula))

    return(validate_LM(res))
  }


  }



# small ball :)
