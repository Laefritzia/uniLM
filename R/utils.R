# here kleine useful functions reingeben

# > parameter checks ----------------------
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

# > MOM helper ----------

# use inside functions with MOM, where conflicts could arise
# when specified an external beta/formula argument, even tho
# these are usually stored inside the LM objects we intend to use for
# these MOM functions.
#
# i hate this function never use it
#
# momOverwrite<-function(var){
#   if(is.null(get(var,envir=parent.frame()))){
#     return(data[[var]])
#   } else{
#     warning(paste0(var, " specified, ignoring value stored in LM-object (data)."))
#     return(get(var, envir=parent.frame()))
#   }
# }



