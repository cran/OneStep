
require(OneStep)
require(actuar)

n <- 1e3

####   G1 <- c("norm", "exp", "lnorm", "invgauss", "pois", "geom") # Distributions for which the explicit MLE is returned ####   

x <- rexp(n)
benchonestep(x, "exp", c("mme", "mle")) 
benchonestep(x, "exp", c("mm", "one")) 


####   G3 <- c("gamma", "beta", "nbinom") # Distributions with a default MME initial guess estimator ####   

x <- rbeta(n, 3, 2)
benchonestep(x, "beta", c("mme", "mle")) 
benchonestep(x, "beta", c("mle", "one")) 
