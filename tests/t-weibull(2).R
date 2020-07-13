

require(OneStep)

n <- 1e3

theta<-c(0.8,3)
o.sample<-rweibull(n,shape = theta[1],scale = 1/theta[2])

#i <- 0
  
dweibull2 <- function(x, shape, scale, log=FALSE)
{
  #cat("\t\tx", i, "\ttheta", shape, "/-/", scale, "\n")
  #print(str(x))
  if(is.function(x))
  {
    
    #print(alist(x))
    #print(formals(x))
    stop("x is a function")
  }
  #i <<- i+1
  dweibull(x = x, shape = shape, scale = scale, log = log)
}
args(dweibull)

onestep(o.sample, "weibull2", method="numeric", start=list(shape=1, scale=1))
mledist(o.sample, "weibull2", start=list(shape=1, scale=1))
