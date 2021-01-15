
require(OneStep)
require(actuar)

n <- 1e3

####   G1 <- c("norm", "exp", "lnorm", "invgauss", "pois", "geom") # Distributions for which the explicit MLE is returned ####   

x <- rexp(n)
os <- onestep(x, "exp", method="numeric")
summary(os)
print(os)
class(os)
plot(os)

x <- rlnorm(n, 0, 1)
os <- onestep(x, "lnorm", method="numeric")
summary(os)
plot(os)

####   G2 <- c("unif", "logis") # Special distributions ####   


x <- rlogis(n)
os <- onestep(x, "logis", method="numeric")
summary(os)
le <- fitdist(x, "logis")
summary(le)


####   G3 <- c("gamma", "beta", "nbinom") # Distributions with a default MME initial guess estimator ####   

x <- rbeta(n, 3, 2)
os <- onestep(x, "beta", method="numeric")
summary(os)
os <- onestep(x, "beta", method="numeric", init=list(shape1=1, shape2=1))
summary(os)

x <- rgamma(n, 3, 2)
os <- onestep(x, "gamma", method="numeric")
summary(os)
os <- onestep(x, "gamma", method="numeric", init=list(shape=1, scale=1.5))
summary(os)
os <- onestep(x, "gamma", method="numeric", init=list(shape=1, rate=1.5))
summary(os)

x <- rnbinom(n, 3, 1/3) #mu=size*(1-prob)/prob=6
os <- onestep(x, "nbinom", method="numeric")
summary(os)
os <- onestep(x, "nbinom", method="numeric", init=list(size=10, mu=2))
summary(os)
os <- onestep(x, "nbinom", method="numeric", init=list(size=10, prob=5/6))
summary(os)




####   G4 <- c("cauchy", "weibull", "pareto") # Distributions with a special initial guess estimator ####   

x <- rweibull(n, 3, 2)
os <- onestep(x, "weibull", method="numeric")
summary(os)
os <- onestep(x, "weibull", method="numeric", init=list(shape=1, scale=1))
summary(os)


x <- rcauchy(n, 1, 2)
os <- onestep(x, "cauchy", method="numeric", control=list(param_t=1/4))
summary(os)
os <- onestep(x, "cauchy", method="numeric", init=list(location=0, scale=1))
summary(os)


x <- rpareto(n, 1, 2)
os <- onestep(x, "pareto", method="numeric")
summary(os)
os <- onestep(x, "pareto", method="numeric", init=list(shape=1, scale=1))
summary(os)
