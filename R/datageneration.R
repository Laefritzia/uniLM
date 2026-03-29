#' generates \(corrupt\) data of class LM for different scenarios
#'
#' @param n integer specifiyng number of observations
#' @param beta true linear relationship between response and predictors for CLEAN data (Inliers)
#' @param O in [0,1]: percentage of data to be corrupted (O..Outliers)
#' @param OO worst range for Outliers
#' @param scenario as specified in seminar-paper
#' @returns named list of class LM
#' @importFrom stats rnorm
#' @export
corrupt_data <- function(n, scenario=c("a", "b", "c"),
                        beta=list(a=c(1,3), b=NULL, c=NULL),
                         O=0.2, OO=100, eps_var=0.3
                        ){

  # To-Do:
  # .) more scenarios: good data, heavy tails, multivariable,...
  # .) deal with intercept

  # Checks
  invisible(
    mapply(check_1num, list(n,O,OO,eps_var),c("num", "prob", rep("pos_num",2)))
  )

  eps<-rnorm(n)*eps_var
  scenario<-match.arg(scenario)

  if(scenario %in% "a"){
  # clean data (Inliers) follows true beta + noise
  # corrupt data (Outliers: cauchy and random): shifted intercept, opposite slope

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

  # building the model
  betaA<-beta$a
  form<-as.formula("y~x2") # formula
  clean<-as.integer(pos%in%"I") # index clean data
  X <- model.matrix(form, data=data.frame(x2=x2,y=numeric(n))) # Intercept+covariates

  y<- clean*(X%*%betaA+eps) + # Clean data
  (1-clean)*(OO+betaA[1]+X[,-1,drop=FALSE]%*%-betaA[-1]+eps) # corrupt data

  res<-new_LM(list(
    data=data.frame(x2=x2,eps=eps,y=y, clean=clean),
    beta=betaA,
    form=form))

  return(validate_LM(res))
  }

  if(scenario %in% "b"){
    return(cat("scenario b to be done"))
  }

  if(scenario %in% "c"){
    return(cat("scenario c to be done"))
  }
}
