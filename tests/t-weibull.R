

require(OneStep)

n <- 1e3
set.seed(1)
theta <- c(0.8,3)
o.sample <- rweibull(n,shape = theta[1], scale = 1/theta[2])

onestep(o.sample, "weibull")
