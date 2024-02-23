
require(OneStep)

#### onestep numerical ####
n <- 1e3
x <- c(rexp(n), rgamma(n, 1/2))
dexpgamma <- function(x, shape1, rate1, rate2)
  1/2*dgamma(x, shape1, rate1) + 1/2*dexp(x, rate2)

onestep(x, "expgamma", control=list(trace=1), start=list(shape1=2, rate1=2, rate2=2))
onestep(x, "expgamma", control=list(trace=2), start=list(shape1=2, rate1=2, rate2=2))
onestep(x, "expgamma", control=list(trace=3), start=list(shape1=2, rate1=2, rate2=2))

onestep(x, "expgamma", control=list(trace=3, REPORT=2), start=list(shape1=2, rate1=2, rate2=2), optim.method="BFGS")


#### onestep closed form with subsample ####

n <- 1e3

x <- 1+2*rt(n, 15)

require(extraDistr)

onestep(x, "lst", control=list(trace=1, delta=1/3))
onestep(x, "lst", control=list(trace=2, delta=1/3))

onestep(x, "lst", control=list(trace=1, delta=2/3))
onestep(x, "lst", control=list(trace=2, delta=2/3))
onestep(x, "lst", control=list(trace=3, delta=2/3))

