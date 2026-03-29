
corrupt_data <- function(n, beta=c(1,3), p=length(beta)-1,
                         O=0.2, OO=100, eps_var=0.3,
                         form=as.formula("y~x2"), covariables=1, ...){
  # O...percentage of data to be corrupted (in sim variieren)
  # OO...worst range for outlier

  # To-Do: auch normale Daten, usw.
  #        deal with intercept

  eps<-rnorm(n)*eps_var

  # drawing (and shuffling) Inliers and Outliers based on O-probabilities
  pos<-factor(
    sample(1:3, n, prob = c(1-O, O/2, O/2), replace=TRUE),
    labels=c("I", "O1", "O2")
  )

  clean<-as.integer(pos%in%"I")

  # generating corrupted vector
  x2<-numeric(n)
  x2[pos%in%"I"]<-rnorm(sum(pos %in% "I")) #normal
  x2[pos%in%"O1"]<-rnorm(sum(pos %in% "O1"))/rnorm(sum(pos %in% "O1")) #cauchy
  x2[pos%in%"O2"]<-sample(seq(-OO,OO,by=0.001), sum(pos %in% "O2"), replace=TRUE) #random

  X <- model.matrix(form, data=data.frame(x2=x2,y=numeric(n)))

  #clean data follows lm
  y<- clean*(X%*%beta+eps) +
  #corrupt data: shifted intercept, opposite slope
  (1-clean)*(OO+beta[1]+X[,-1,drop=FALSE]%*%-beta[-1]+eps)

  # und nicht hardcoden, sondern dynamisch die variablen rein
  test<-data.frame(x2=x2,eps=eps,y=y, clean=clean)

  return(list(data=test,beta=beta))
}
