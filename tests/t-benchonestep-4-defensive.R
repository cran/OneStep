
require(OneStep)
require(actuar)

n <- 1e3

x <- rexp(n)
try( benchonestep(x, "exp", "jaimeleraisindetable") )
try( benchonestep(x, "exp") )
try( benchonestep(x, "exp2", "mle") )
try( benchonestep(matrix(x, 2, n/2), "exp2", "mle") )


