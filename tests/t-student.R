
require(OneStep)

n <- 1e3

#### standard student distribution ####

theta <- c("df"=3/2, "ncp"=0)
o.sample <- rt(n, df = theta["df"], ncp = theta["ncp"])

mledist(o.sample, "t", start=list("df"=2), fix.arg=list("ncp"=0))
onestep(o.sample, "t", control=list(delta=0.3))
onestep(o.sample, "t", control=list(delta=0.75))

