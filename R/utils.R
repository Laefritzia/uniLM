# here kleine useful functions reingeben

# useful for parameter checks
check_1num <-function(x, type=c("int", "intExact", "num", "prob", "pos_num")){
  t<-match.arg(type)
  is.good<-switch(t,
                  "int"=function(i) i%%1==0,
                  "intExact"=is.integer,
                  "num"=is.numeric,
                  "prob"=function(i) i>=0&i<=1,
                  "pos_num"=function(i) i>=0)

  if(!is.good(x)||length(x)!=1){
    stop(paste0(deparse(substitute(x)), " must be a single number of type ", t))
  }
}


