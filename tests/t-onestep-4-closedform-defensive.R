
require(OneStep)

n <- 1e3


x <- rexp(n)
try(onestep(matrix(x, 4, n/4), "exp"))

try(onestep(x[1], "exp"))

try(onestep(x, "exp", method="jaimeleraisindetable"))

try(onestep(x, "jaimeleraisindetable"))

try(onestep(x, "exp", init=list(jaimeleraisindetable=1)))

try(onestep(x, "exp", init=list(rate=1)))    

try(onestep(x, "gamma", init=list(rate=1)))    

try(onestep(x, "gamma", init=list(shape=1)))    
