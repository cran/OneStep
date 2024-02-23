

require(OneStep)

n <- 1e3

theta <- c(shape=0.8, scale=3)
o.sample <- rweibull(n, shape=theta["shape"], scale=theta["scale"])

#### numerical one-step with sub-sample init ####
#i <- 0
  
dweibull2 <- function(x, shape, scale, log=FALSE)
{
  #cat("\t\tx", i, "\ttheta", shape, "/-/", scale, "\n")
  if(is.function(x))
  {
    stop("x is a function")
  }
  #i <<- i+1
  dweibull(x = x, shape = shape, scale = scale, log = log)
}
args(dweibull)

onestep(o.sample, "weibull2", method="numeric", start=list(shape=1, scale=1), control=list(trace=1))

mledist(o.sample, "weibull2", start=list(shape=1, scale=1))

#### numerical one-step with user-supplied ####

onestep(o.sample, "weibull", method="numeric", start=list(shape=1, scale=1), control=list(trace=1))
