
require(OneStep)

n <- 1e3

#### standard student distribution ####

#theta <- c("df"=3/2, "ncp"=1)
#o.sample <- rt(n, df = theta["df"], ncp = theta["ncp"])

#ML1 <- mledist(o.sample, "t", start=list("df"=2, "ncp"=0))
#onestep(o.sample, "t", control=list(delta=1/2))
#onestep(o.sample, "t", control=list(delta=.75))

#### location-scaled student distribution ####

require(extraDistr)

curve(dlst(x, 3/2), -10, 10)

theta <- c("df"=3/2, "mu"=1, "sigma"=2)
o.sample <- rlst(n, df=theta["df"], mu=theta["mu"], sigma=theta["sigma"])

ML1 <- fitdist(o.sample, "lst", start=list("df"=10, "mu"=median(o.sample), "sigma"=IQR(o.sample)/2))
OS1 <- onestep(o.sample, "lst", control=list(delta=0.4))
OS2 <- onestep(o.sample, "lst", control=list(delta=0.8))


cbind("theo"=theta, "MLE"=coef(ML1)[names(theta)], 
      "OS1"=coef(OS1)[c("df", "mu", "sigma")],"OS2"=coef(OS2)[c("df", "mu", "sigma")])
