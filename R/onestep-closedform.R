onestep_closedformula  <-  function(data, distname, control, ...)
{
  
  n  <-  length(data)
  delta <- control$delta
  
  if (distname == "norm") { #MLE=ME already optimal
    m  <-  mean(data)
    v  <-  (n - 1)/n * var(data)
    estimate  <-  c(mean = m, sd = sqrt(v))
    order  <-  1:2
  }
  else if (distname == "exp") { #MLE=ME already optimal
    m  <-  mean(data)
    estimate  <-  c(rate = 1/m)
    order  <-  1
  }
  else if (distname == "lnorm") { #MLE already  optimal 
    if (any(data <= 0)) 
      stop("values must be positive to fit a lognormal distribution")
    meanlog  <-  mean(log(data))
    sdlog    <-  sqrt(mean((log(data)-meanlog)^2))
    estimate  <-  c(meanlog = meanlog, sdlog = sdlog)
    order  <-  0  #ME non optimal !
  }
  else if(distname=="invgauss"){ #MLE already  optimal 
    m  <-  mean(data)
    d  <-  mean(1/data)-1/m
    estimate  <-  c(mean = m, dispersion = d)
    order  <-  0  #ME non optimal !
  }
  else if (distname == "pois") { #MLE=ME  already optimal
    m  <-  mean(data)
    estimate  <-  c(lambda = m)
    order  <-  1
  }
  else if (distname == "geom") { #MLE=ME already  optimal
    m  <-  mean(data)
    prob  <-  if (m > 0) 
      1/(1 + m)
    else NaN
    estimate  <-  c(prob = prob)
    order  <-  1
  }
  else if (distname == "gamma") { #Done
    m  <-  mean(data)
    v  <-  (n - 1)/n * var(data)
    shape  <-  m^2/v
    rate  <-  m/v
    
    I <- matrix(0,2,2)
    I[1,1] <- trigamma(shape)
    I[2,1] <- -1/rate
    I[1,2] <- -1/rate
    I[2,2] <- shape/rate^2
    
    Score <- matrix(c(sum(log(rate)-digamma(shape)+log(data)),sum(shape/rate-data)),
                  2,
                  1)
    LCE <- c(shape,rate)+ 1/n*solve(I,Score)
    
    estimate  <-  c(shape = LCE[1], rate = LCE[2])
    order  <-  1:2
  }
  else if (distname == "beta") { #Done
    if (any(data < 0) | any(data > 1)) 
      stop("values must be in [0-1] to fit a beta distribution")
    m  <-  mean(data)
    v  <-  (n - 1)/n * var(data)
    aux  <-  m * (1 - m)/v - 1
    shape1  <-  m * aux
    shape2  <-  (1 - m) * aux
    
    I <- matrix(0,2,2)
    tritmp <- trigamma(shape1+shape2)
    I[1,1] <- trigamma(shape1)-tritmp
    I[2,1] <- -tritmp
    I[1,2] <- -tritmp
    I[2,2] <- trigamma(shape2)-tritmp
    
    ditmp <- digamma(shape1+shape2)
    Score <- matrix(c(sum(ditmp-digamma(shape1)+log(data)),sum(ditmp-digamma(shape2)+log(1-data))),2,1)
    LCE <- c(shape1,shape2)+ 1/n*solve(I,Score)
    
    estimate  <-  c(shape1 = LCE[1], shape2 = LCE[2])
    order  <-  1:2
    
  }
  else if (distname == "cauchy"){ #To finish
    
    param_t <- control$param_t 
    if(is.null(param_t))
      param_t <- 1/4
    
    Phi <- mean(exp(1i*param_t*data))
    S <- Re(Phi)
    Z <- Im(Phi)
    m <- atan(Z/S)/param_t
    d <- -log(sqrt(S^2+Z^2))/param_t
    
    denom <- d^2+(data-m)^2
    LCE <- c(m,d) + 2*d^2 *c(mean(2*(data-m)/denom),mean(1/d-2*d/denom))
    
    estimate  <-  c(location = LCE[1], scale = LCE[2])
    order  <- 0
    
  }
  else if (distname == "weibull") {
    Score <- matrix(0,2,1)
    I <- matrix(0,2,2)
    if(is.null(delta))
      stop("wrong delta parameter")
    ndelta <- as.integer(n^(delta))
    idxsubsample <- sample(1:n, ndelta, replace=FALSE)
    
    esttmp  <-  mledist(data[idxsubsample], distr="weibull",...)$estimate
    
    shape  <-  esttmp[1]
    scale  <-  esttmp[2]
    rate   <-  1/scale
    
    shape  <-  esttmp[1]
    scale  <-  esttmp[2]
    rate   <-  1/scale
    
    Score[1,1] <- mean(1/shape+ log(rate)+log(data)-log(rate*data)*(rate*data)^shape)
    Score[2,1] <- mean(shape/rate-shape*(data^shape)*rate^(shape-1))
    
    I[1,1] <-  1/shape^2*(trigamma(1)+digamma(2))
    I[2,1] <-  1/rate*digamma(2) 
    I[1,2] <-  1/rate*digamma(2)  
    I[2,2] <-  shape^2/rate^2
    
    LCE  <-  c(shape,rate) + solve(I,Score)  
    estimate  <-  c(shape = LCE[1], scale = 1/LCE[2])
    order  <-  0
  } 
  else if (distname == "pareto") {
    if(is.null(delta))
      stop("wrong delta parameter")
    ndelta <- as.integer(n^(delta))
    idxsubsample <- sample(1:n, ndelta, replace=FALSE)
    
    esttmp<-mledist(data[idxsubsample],"pareto",...)$estimate
    
    shape <- esttmp[1]
    scale <- esttmp[2]
    
    Score<- matrix(0,2,1)
    Score[1,1]<- mean(1/shape+log(scale)-log(scale+data))
    Score[2,1]<- mean(shape/scale-(shape+1)*(data+scale)^(-1))
    
    I<-matrix(0,2,2)
    I[1,1]<- 1/(shape^2)
    I[2,1]<- -1/(scale*(shape+1))
    I[1,2]<- -1/(scale*(shape+1)) 
    I[2,2]<- shape/(scale^2*(shape+2))
    
    LCE<-c(shape,scale)+ solve(I,Score)
    estimate  <-  c(shape = LCE[1], scale = LCE[2])
    order  <-  0
  }  
  else if (distname == "nbinom") {
    m  <-  mean(data)
    v  <-  (n - 1)/n * var(data)
    size  <-  if (v > m) 
      r <- m^2/(v - m)
    else NaN
    mu <- m
    
    Ihat <- matrix(0,2,2)
    Ihat[1,1] <- -mean(trigamma(r+data)-trigamma(r)+1/r-1/(mu+r)+(data-mu)/(mu+r)^2)
    Ihat[1,2] <- -mean(r/(mu+r)^2-1/(mu+r)+data/(mu+r)^2)
    Ihat[2,1] <- Ihat[1,2]
    Ihat[2,2] <- -mean(r/(mu+r)^2+data*(1/(mu+r)^2-1/mu^2))
    
    sc1 <- mean(digamma(r+data)-digamma(r)+1+log(r)-r/(mu+r)-log(mu+r)-data/(mu+r))
    sc2 <- mean(-r/(mu+r)+data*(1/mu-1/(mu+r)))
    
    LCE <- c(r,mu)+ solve(Ihat,matrix(c(sc1,sc2),2,1))
    
    estimate  <-  c(size = LCE[1], mu = LCE[2])
    order  <-  1:2
  }
  else if (distname == "unif") {
    m  <-  mean(data)
    v  <-  (n - 1)/n * var(data)
    min1  <-  m - sqrt(3 * v)
    max1  <-  m + sqrt(3 * v)
    estimate  <-  c(min1, max1)
    order  <-  1:2
  }
  else if (distname == "logis") {
    m  <-  mean(data)
    v  <-  (n - 1)/n * var(data)
    scale  <-  sqrt(3 * v)/pi
    estimate  <-  c(location = m, scale = scale)
    order  <-  1:2
  }else
  {
    stop("unsupported distribution")
  }
  
  #if (exists(ddistname)) 
  #    loglikval  <-  loglik(estimate, fix.arg, data, ddistname)
  #else loglikval  <-  NULL
  loglikval <- 0
  
  list(estimate = estimate, convergence = 0, value = NULL, 
              hessian = NULL, optim.function = NULL, opt.meth = NULL, 
              fix.arg = NULL, fix.arg.fun = NULL, weights = NULL, 
              counts = NULL, optim.message = NULL, loglik = loglikval, 
              method = "closed formula", order = order, memp = NULL)
}