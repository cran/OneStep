onestep_closedformula <- function(data, distname, control, init,...)
{
  
  n <- length(data)
  
  G1 <- c("norm", "exp", "lnorm", "invgauss", "pois", "geom") # Distributions for which the explicit MLE is returned
  G2 <- c("unif", "logis") # Special distributions
  G3 <- c("gamma", "beta", "nbinom","chisq") # Distributions with a default MME initial guess estimator
  G4 <- c("cauchy", "weibull", "pareto") # Distributions with a special initial guess estimator
  
  if (!is.element(distname, c(G1,G2,G3,G4)))
  {
    stop("unsupported distribution")
  }
  
  if(!missing(init)){
    if (is.element(distname, c(G1,G2))){ 
      warning("The initial estimate is not taken into account for this distribution")
      initb <- FALSE
    }else{
      initb <- TRUE
      #test sur le type et la dimension de init
    }
  }else{
    initb <- FALSE
  }
  
  ### G1 - Distributions for which the explicit MLE is returned
  
  if (distname == "norm") { #MLE=MME already optimal
    m <- mean(data)
    v <- (n - 1)/n * var(data)
    estimate <- c(mean = m, sd = sqrt(v))
    order <- 1:2
  }
  
  if (distname == "exp") { #MLE=MME already optimal
    m <- mean(data)
    estimate <- c(rate = 1/m)
    order <- 1
  }
  
  if (distname == "lnorm") { #MLE already  optimal 
    if (any(data <= 0)) 
      stop("values must be positive to fit a lognormal distribution")
    meanlog <- mean(log(data))
    sdlog   <- sqrt(mean((log(data)-meanlog)^2))
    estimate <- c(meanlog = meanlog, sdlog = sdlog)
    order <- 0  #ME non optimal !
  }
  
  if(distname=="invgauss"){ #MLE already  optimal 
    m <- mean(data)
    d <- mean(1/data)-1/m
    estimate <- c(mean = m, dispersion = d)
    order <- 0  #ME non optimal !
  }
  
  if (distname == "pois") { #MLE=MME  already optimal
    m <- mean(data)
    estimate <- c(lambda = m)
    order <- 1
  }
  
  if (distname == "geom") { #MLE=MME already  optimal
    m <- mean(data)
    prob <- if (m > 0) 
      1/(1 + m)
    else NaN
    estimate <- c(prob = prob)
    order <- 1
  }
  
  ### G2 - Special Distributions
  
  if (distname == "unif") { #should be in G1
    m <- mean(data)
    v <- (n - 1)/n * var(data)
    min1 <- m - sqrt(3 * v)
    max1 <- m + sqrt(3 * v)
    estimate <- c(min1, max1)
    order <- 1:2
  }
  
  if (distname == "logis") { #should be in G3, see also section 3.3 of Handbook of logistic distr, p64
    m <- mean(data)
    v <- (n - 1)/n * var(data)
    scale <- sqrt(3 * v)/pi
    estimate <- c(location = m, scale = scale)
    order <- 1:2
  }
  
  ### G3 - Distributions with a default MME initial guess estimator 
  
  if (distname == "gamma") { #Done
    
    if (initb){
      shape <- init$shape
      rate <- ifelse(is.null(init$rate), 1/init$scale, init$rate)
    }else{
      m <- mean(data)
      v <- (n - 1)/n * var(data)
      shape <- m^2/v
      rate <- m/v
    }
    
    IFisher <- matrix(0,2,2)
    IFisher[1,1] <- trigamma(shape)
    IFisher[2,1] <- -1/rate
    IFisher[1,2] <- -1/rate
    IFisher[2,2] <- shape/rate^2
    
    Score <- matrix(c(sum(log(rate)-digamma(shape)+log(data)),
                      sum(shape/rate-data)), 2, 1)
    LCE <- c(shape,rate)+ 1/n*solve(IFisher,Score)
    
    estimate <- c(shape = LCE[1], rate = LCE[2])
    order <- 1:2
  }
  
  if (distname == "beta") { #Done
    if (any(data < 0) | any(data > 1)) 
      stop("values must be in [0-1] to fit a beta distribution")
    
    
    if(initb){
      shape1 <- init$shape1
      shape2 <- init$shape2
    }else{
      m <- mean(data)
      v <- (n - 1)/n * var(data)
      aux <- m * (1 - m)/v - 1
      shape1 <- m * aux
      shape2 <- (1 - m) * aux
    }
    
    IFisher <- matrix(0,2,2)
    tritmp <- trigamma(shape1+shape2)
    IFisher[1,1] <- trigamma(shape1)-tritmp
    IFisher[2,1] <- -tritmp
    IFisher[1,2] <- -tritmp
    IFisher[2,2] <- trigamma(shape2)-tritmp
    
    ditmp <- digamma(shape1+shape2)
    Score <- matrix(c(sum(ditmp-digamma(shape1)+log(data)),sum(ditmp-digamma(shape2)+log(1-data))),2,1)
    LCE <- c(shape1,shape2)+ 1/n*solve(IFisher,Score)
    
    estimate <- c(shape1 = LCE[1], shape2 = LCE[2])
    order <- 1:2
    
  }
  
  if (distname == "nbinom") {
    
    if (initb){
      size <- init$size
      mu <- ifelse(is.null(init$mu), size/init$prob-size, init$mu)
    }else{
      m <- mean(data)
      v <- (n - 1)/n * var(data)
      r <- m^2/(v - m)
      size <- ifelse(v > m, r, NaN)
      mu <- m
    }
    musize1 <- 1/(mu+size)
    musize2 <- 1/(mu+size)^2
    
    Ihat <- matrix(0,2,2)
    Ihat[1,1] <- -mean(trigamma(size+data)-trigamma(size) + 1/size - musize1 + (data-mu)*musize2)
    Ihat[1,2] <- -mean(size*musize2 - 1/(mu+size) + data*musize2)
    Ihat[2,1] <- Ihat[1,2]
    Ihat[2,2] <- -mean(size*musize2 + data*(musize2 - 1/mu^2))
    
    sc1 <- mean(digamma(size+data)-digamma(size) + 1 + log(size) - (size+data)*musize1-log(mu+size)) 
    sc2 <- mean(-size*musize1 + data*(1/mu-musize1))
    
    LCE <- c(size,mu) + solve(Ihat,matrix(c(sc1,sc2),2,1))
    
    estimate <- c(size = LCE[1], mu = LCE[2])
    order <- 1:2
  }
  
  if (distname=="chisq"){
    
    if (initb){
      df <- init$df
    }else{
      df<-mean(data)
    }
    
    Ihat<-trigamma(df/2)/4
    sc<-mean(log(data/2)/2-digamma(df/2)/2)
    
    LCE<-df+ sc/Ihat

    estimate <- c(df=LCE)
    order <- 1
  }
  
  ### G4 - Distributions with a special initial guess estimator
  
  delta <- control$delta
  
  if (distname == "cauchy"){
    
    if(initb){
      mylocation <- init$location
      myscale <- init$scale
    }else{
      #Quantiles Matching
      mylocation <- as.numeric(quantile(data,1/2,type=3))
      myscale <- as.numeric(diff(quantile(data,c(1/4,3/4),type=3)))/2  
    }
  
    denom <- myscale^2+(data-mylocation)^2
    LCE <- c(mylocation, myscale) + 2*myscale^2 *c(mean(2*(data-mylocation)/denom),
                                                   mean(1/myscale-2*myscale/denom))
    
    estimate <- c(location = LCE[1], scale = LCE[2])
    order  <- 0
    
  }
  
  if (distname == "weibull") {
   
    if (initb){
      
      shape<-init$shape
      scale<-init$scale
      #dweibull in stats pkg does not allow for a rate parameter
      rate<-1/scale 
      
    }else{
      
      #if(is.null(delta))
      #  stop("wrong delta parameter")
      #ndelta <- as.integer(n^(delta))
      #idxsubsample <- sample(1:n, ndelta, replace=FALSE)
      #esttmp <- mledist(data[idxsubsample], distr="weibull",...)$estimate
    
      #shape <- esttmp[1]
      #scale <- esttmp[2]
      #dweibull in stats pkg does not allow for a rate parameter

      #rate  <- 1/scale
      
      Y<- log(sort(data))
      X<- X<- log(-log(1-(1:n)/(n+1))) #voir aussi ASTM Procedure
      res<-as.numeric(lm(Y~X)$coef)
      shape<-1/res[2]
      rate<-1/exp(res[1])
      
    }
    
    Score <- matrix(0,2,1)
    Score[1,1] <- mean(1/shape+ log(rate)+log(data)-log(rate*data)*(rate*data)^shape)
    Score[2,1] <- mean(shape/rate-shape*(data^shape)*rate^(shape-1))
    
    IFisher <- matrix(0,2,2)
    IFisher[1,1] <-  1/shape^2*(trigamma(1)+digamma(2)^2)
    IFisher[2,1] <-  1/rate*digamma(2) 
    IFisher[1,2] <-  1/rate*digamma(2)  
    IFisher[2,2] <-  shape^2/rate^2
    
    LCE <- c(shape,rate) + solve(IFisher,Score) 
    #dweibull in stats pkg does not allow for a rate parameter
    estimate <- c(shape = LCE[1], scale = 1/LCE[2])
    order <- 0
  } 
  
  if (distname == "pareto") {
    
    if (initb){
      shape<-init$shape
      scale<-init$scale
    }else{
      
      if(is.null(delta))
        stop("wrong delta parameter")
      ndelta <- as.integer(n^(delta))
      idxsubsample <- sample(1:n, ndelta, replace=FALSE)
    
      esttmp <- mledist(data[idxsubsample],"pareto",...)$estimate
    
      shape <- esttmp[1]
      scale <- esttmp[2]
    }
    
    Score<- matrix(0,2,1)
    Score[1,1] <-  mean(1/shape+log(scale)-log(scale+data))
    Score[2,1] <-  mean(shape/scale-(shape+1)*(data+scale)^(-1))
    
    IFisher <- matrix(0,2,2)
    IFisher[1,1] <- 1/(shape^2)
    IFisher[2,1] <- -1/(scale*(shape+1))
    IFisher[1,2] <- -1/(scale*(shape+1)) 
    IFisher[2,2] <- shape/(scale^2*(shape+2))
    
    LCE<-c(shape,scale)+ solve(IFisher,Score)
    estimate <- c(shape = LCE[1], scale = LCE[2])
    order <- 0
  }  
  
  list(estimate = estimate, convergence = 0, value = NULL, 
              hessian = NULL, optim.function = NULL, opt.meth = NULL, 
              fix.arg = NULL, fix.arg.fun = NULL, weights = NULL, 
              counts = NULL, optim.message = NULL, loglik = NULL, 
              method = "closed formula", order = order, memp = NULL)
}