
# > 1 ------------------
new_LM <- function(res){ # can pass args, eg strat_var
  obj <- structure(res, class="LM")
  #attributes(obj, "strat_var") <- strat_var
}

# > 2 ------------------
LM <- function(formula, data){

  # To-DO:
  # kommentieren
  # checks und eingabe
  # sollten mit intercept dealen, hier Annahme dass intercept dabei
  # maybe auch andere robuste estimators(m) und varianzen rein?s

  p<-length(attr(terms(formula),"term.labels")) # p sind predictoren ohne intercept, wie bei Krivo Elmod VO 1
  n<-nrow(data)

  X<-model.matrix(formula, data=data)
  H<-X%*%solve(t(X)%*%X)%*%t(X)   # Projektionsmatrix auf span(X)
  M<- diag(n)-H                   # orthogonales Komplement zu H
  ols<-solve(t(X)%*%X)%*%t(X)%*%y # ordinary-least-square-estimator

  yhat<-X%*%ols # H%*%y
  ehat<-y-yhat  # residuals

  # standardized residuals
  shat<-sum(ehat^2)/(n-(p+1))
  r_std<-ehat/(sqrt(shat*(1-diag(H))))
  # studentized residuals (umwandlung: https://online.stat.psu.edu/stat462/node/247/)
  r_t<- r_std*sqrt((n-p-2)/(n-p-1-r_std^2))
  # Cooks DIstance (https://de.wikipedia.org/wiki/Cook-Abstand)
  cook<-r_t^2/(p+1)*(diag(H)/(1-diag(H)))

  return(new_LM(
    list()# was auch immer wir als output reingeben wollen
  ))


}

# > 3 ------------------
residuals.LM <- function(){}

# > 4 ------------------
plot.LM <- function(){
  # vielleicht nur ein wrapper, der auf diagnostics zugreifft?
}


# > 5 ------------------
confint.LM <- function(object, data, parm = "trt", level = 0.95,
                           test=c("Score", "Wald"),
                           tol=10^-9, maxiter=20, ...){

      if (!inherits(object, "my_cox")){
        stop("Object must be class 'my_cox'.")
      }

      test <- match.arg(test)
    }






