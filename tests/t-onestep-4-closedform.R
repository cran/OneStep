
require(OneStep)
require(actuar)

n <- 1e3

####   G1 <- c("norm", "exp", "lnorm", "invgauss", "pois", "geom") # Distributions for which the explicit MLE is returned ####   

x <- rexp(n)
os <- onestep(x, "exp")
summary(os)
print(os)
class(os)
plot(os)

x <- rlnorm(n, 0, 1)
os <- onestep(x, "lnorm", method="closed")
summary(os)
plot(os)

####   G2 <- c("unif", "logis") # Special distributions ####   

x <- runif(n)
os <- onestep(x, "unif")
summary(os)
le <- fitdist(x, "unif")
summary(le)


x <- rlogis(n)
os <- onestep(x, "logis")
summary(os)
le <- fitdist(x, "logis")
summary(le)


####   G3 <- c("gamma", "beta", "nbinom") # Distributions with a default MME initial guess estimator ####   

x <- rbeta(n, 3, 2)
os <- onestep(x, "beta")
summary(os)
os <- onestep(x, "beta", init=list(shape1=1, shape2=1))

x <- rgamma(n, 3, 2)
os <- onestep(x, "gamma", init=list(shape=1, scale=1))
summary(os)
os <- onestep(x, "gamma", init=list(shape=10, rate=1))
summary(os)

x <- rnbinom(n, 3, 1/3) #mu=size*(1-prob)/prob=6
os <- onestep(x, "nbinom")
summary(os)
os <- onestep(x, "nbinom", init=list(size=10, mu=2))
summary(os)
os <- onestep(x, "nbinom", init=list(size=10, prob=1/2))
summary(os)




####   G4 <- c("cauchy", "weibull", "pareto") # Distributions with a special initial guess estimator ####   

x <- rweibull(n, 3, 2)
os <- onestep(x, "weibull")
summary(os)
os <- onestep(x, "weibull", init=list(shape=1, scale=1))
summary(os)


x <- rcauchy(n, 1, 2)
os <- onestep(x, "cauchy", control=list(param_t=1/4))
summary(os)
os <- onestep(x, "cauchy", init=list(location=0, scale=1))
summary(os)


x <- rpareto(n, 1, 2)
os <- onestep(x, "pareto")
summary(os)
os <- onestep(x, "pareto", init=list(shape=1, scale=1))
summary(os)
