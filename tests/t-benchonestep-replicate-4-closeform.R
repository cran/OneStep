
require(OneStep)


n <- 1e3

####   G1 <- c("norm", "exp", "lnorm", "invgauss", "pois", "geom") # Distributions for which the explicit MLE is returned ####   


benchonestep.replicate(n, 4, "exp", rate=1)
benchonestep.replicate(n, 4, "exp", rate=1, methods="mle")
  




####   G3 <- c("gamma", "beta", "nbinom") # Distributions with a default MME initial guess estimator ####   

benchonestep.replicate(n, 4, "beta", shape1=1/2, shape2=3/2)
if(FALSE)
  benchonestep.replicate(n, 4, "beta", shape1=1/2, shape2=3/2, ncpus=2)