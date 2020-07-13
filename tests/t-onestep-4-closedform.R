
require(OneStep)

n <- 1e3


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

x <- rcauchy(n)
os <- onestep(x, "cauchy", control=list(param_t=1/4))

