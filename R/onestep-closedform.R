onestep_closedformula <- function(obs, distname, control, init,...)
{
  n <- length(obs)
  #group fo distributions
  G1 <- c("norm", "exp", "lnorm", "invgauss", "pois", "geom") # Distributions for which the explicit MLE is returned
  G2 <- c("unif", "logis") # Special distributions
  G3 <- c("gamma", "beta", "nbinom","chisq") # Distributions with a default MME initial guess estimator
  G4 <- c("cauchy", "weibull", "pareto", "t", "lst") # Distributions with a special initial guess estimator
  
  #some optim control
  con <- list(trace = 0, fnscale = 1, maxit = 100L, abstol = -Inf, 
              reltol = sqrt(.Machine$double.eps), alpha = 1, beta = 0.5, 
              gamma = 2, REPORT = 10, 
              warn.1d.NelderMead = FALSE, #default value is TRUE in optim()
              type = 1, lmm = 5, factr = 1e+07, pgtol = 0, tmax = 10, temp = 10)
  con[names(control)] <- control
  con$trace <- max(con$trace-1, 0) #decrease trace for optim() only
  
  #number of step of Newton method : no default 
  nbstep <- NULL
  
  if (!is.element(distname, c(G1, G2, G3, G4)))
  {
    stop("unsupported distribution")
  }
  
  if(!missing(init))
  {
    if (is.element(distname, c(G1,G2)))
    { 
      warning("The initial estimate is not taken into account for this distribution")
      initb <- FALSE
    }else
    {
      initb <- TRUE
      #test sur le type et la dimension de init
    }
  }else
  {
    initb <- FALSE
  }
  
  ### Group 1 - Distributions for which the explicit MLE is returned
  
  if (distname == "norm") 
  { #MLE=MME already optimal
    m <- mean(obs)
    v <- (n - 1)/n * var(obs)
    estimate <- c(mean = m, sd = sqrt(v))
    order <- 1:2
    loglik <- NULL
    nbstep <- 0
  }
  
  if (distname == "exp") 
  { #MLE=MME already optimal
    m <- mean(obs)
    estimate <- c(rate = 1/m)
    order <- 1
    loglik <- NULL
    nbstep <- 0
  }
  
  if (distname == "lnorm") 
  { #MLE already  optimal 
    if (any(obs <= 0)) 
      stop("values must be positive to fit a lognormal distribution")
    meanlog <- mean(log(obs))
    sdlog   <- sqrt(mean((log(obs)-meanlog)^2))
    estimate <- c(meanlog = meanlog, sdlog = sdlog)
    order <- 0  #ME non optimal !
    loglik <- NULL
    nbstep <- 0
  }
  
  if(distname=="invgauss")
  { #MLE already  optimal 
    m <- mean(obs)
    d <- mean(1/obs)-1/m
    estimate <- c(mean = m, dispersion = d)
    order <- 0  #ME non optimal !
    loglik <- NULL
    nbstep <- 0
  }
  
  if (distname == "pois") 
  { #MLE=MME  already optimal
    m <- mean(obs)
    estimate <- c(lambda = m)
    order <- 1
    loglik <- NULL
    nbstep <- 0
  }
  
  if (distname == "geom") 
  { #MLE=MME already  optimal
    m <- mean(obs)
    prob <- if (m > 0) 
      1/(1 + m)
    else NaN
    estimate <- c(prob = prob)
    order <- 1
    loglik <- NULL
    nbstep <- 0
  }
  
  ### Group 2 - Special Distributions
  
  if (distname == "unif") 
  { #should be in G1
    m <- mean(obs)
    v <- (n - 1)/n * var(obs)
    min1 <- m - sqrt(3 * v)
    max1 <- m + sqrt(3 * v)
    estimate <- c(min1, max1)
    order <- 1:2
    loglik <- NULL
    nbstep <- 0
  }
  
  if (distname == "logis") 
  { #should be in G3, see also section 3.3 of Handbook of logistic distr, p64
    m <- mean(obs)
    v <- (n - 1)/n * var(obs)
    scale <- sqrt(3 * v)/pi
    estimate <- c(location = m, scale = scale)
    order <- 1:2
    loglik <- NULL
    nbstep <- 0
  }
  
  ### Group 3 - Distributions with a default MME initial guess estimator 
  
  if (distname == "gamma") 
  { #Done
    if (initb)
    {
      shape <- init$shape
      rate <- ifelse(is.null(init$rate), 1/init$scale, init$rate)
      if(control$trace > 1)
      {
        cat("**\t user-supplied init\n")
        print(c(shape, rate))
      }
    }else
    {
      m <- mean(obs)
      v <- (n - 1)/n * var(obs)
      shape <- m^2/v
      rate <- m/v
      if(control$trace > 1)
      {
        cat("**\t default init\n")
        print(c(shape, rate))
      }
    }
    IFisher <- matrix(0,2,2)
    IFisher[1,1] <- trigamma(shape)
    IFisher[2,1] <- -1/rate
    IFisher[1,2] <- -1/rate
    IFisher[2,2] <- shape/rate^2
    
    Score <- matrix(c(sum(log(rate)-digamma(shape)+log(obs)),
                      sum(shape/rate-obs)), 2, 1)
    LCE <- c(shape, rate)+ 1/n*solve(IFisher, Score)
    
    estimate <- c(shape = LCE[1], rate = LCE[2])
    order <- 1:2
    loglik <- NULL
    nbstep <- 1
  }
  
  if (distname == "beta") 
  { #Done
    if (any(obs < 0) | any(obs > 1)) 
      stop("values must be in [0-1] to fit a beta distribution")
    if(initb)
    {
      shape1 <- init$shape1
      shape2 <- init$shape2
      if(control$trace > 1)
      {
        cat("**\t user supplied init\n")
        print(c(shape1, shape2))
      }
    }else
    {
      m <- mean(obs)
      v <- (n - 1)/n * var(obs)
      aux <- m * (1 - m)/v - 1
      shape1 <- m * aux
      shape2 <- (1 - m) * aux
      if(control$trace > 1)
      {
        cat("**\t default init\n")
        print(c(shape1, shape2))
      }
    }
    IFisher <- matrix(0, 2, 2)
    tritmp <- trigamma(shape1+shape2)
    IFisher[1,1] <- trigamma(shape1)-tritmp
    IFisher[2,1] <- -tritmp
    IFisher[1,2] <- -tritmp
    IFisher[2,2] <- trigamma(shape2)-tritmp
    
    ditmp <- digamma(shape1+shape2)
    Score <- matrix(c(sum(ditmp-digamma(shape1)+log(obs)), 
                      sum(ditmp-digamma(shape2)+log(1-obs))), 2, 1)
    LCE <- c(shape1, shape2) + 1/n*solve(IFisher, Score)
    
    estimate <- c(shape1 = LCE[1], shape2 = LCE[2])
    order <- 1:2
    loglik <- NULL
    nbstep <- 1
  }
  
  if (distname == "nbinom") 
  {
    if (initb)
    {
      size <- init$size
      mu <- ifelse(is.null(init$mu), size/init$prob-size, init$mu)
      if(control$trace > 1)
      {
        cat("**\t user supplied init\n")
        print(c(size, mu))
      }
    }else
    {
      m <- mean(obs)
      v <- (n - 1)/n * var(obs)
      r <- m^2/(v - m)
      size <- ifelse(v > m, r, NaN)
      mu <- m
      if(control$trace > 1)
      {
        cat("**\t default init\n")
        print(c(size, mu))
      }
    }
    musize1 <- 1/(mu+size)
    musize2 <- 1/(mu+size)^2
    
    Ihat <- matrix(0,2,2)
    Ihat[1,1] <- -mean(trigamma(size+obs)-trigamma(size) + 1/size - musize1 + (obs-mu)*musize2)
    Ihat[1,2] <- -mean(size*musize2 - 1/(mu+size) + obs*musize2)
    Ihat[2,1] <- Ihat[1,2]
    Ihat[2,2] <- -mean(size*musize2 + obs*(musize2 - 1/mu^2))
    
    sc1 <- mean(digamma(size+obs)-digamma(size) + 1 + log(size) - (size+obs)*musize1-log(mu+size)) 
    sc2 <- mean(-size*musize1 + obs*(1/mu-musize1))
    Score <- matrix(c(sc1,sc2), 2, 1)
    
    LCE <- c(size,mu) + solve(Ihat, Score)
    
    estimate <- c(size = LCE[1], mu = LCE[2])
    order <- 1:2
    loglik <- NULL
    nbstep <- 1
  }
  
  if (distname=="chisq")
  {
    if (initb)
    {
      df <- init$df
      if(control$trace > 1)
      {
        cat("**\t user supplied init\n")
        print(df)
      }
    }else
    {
      df <- mean(obs)
      if(control$trace > 1)
      {
        cat("**\t default init\n")
        print(df)
      }
    }
    
    Ihat <- trigamma(df/2)/4
    sc <- mean(log(obs/2)/2 - digamma(df/2)/2)
    
    LCE <- df + sc/Ihat
    
    estimate <- c(df=LCE)
    order <- 1
    loglik <- NULL
    nbstep <- 1
  }
  
  ### Group 4 - Distributions with a special initial guess estimator
  
  delta <- control$delta
  if (distname == "cauchy")
  {
    if(initb)
    {
      mylocation <- init$location
      myscale <- init$scale
      if(control$trace > 1)
      {
        cat("**\t user supplied init\n")
        print(c(scale=as.numeric(myscale), location=as.numeric(mylocation)))
      }
    }else
    {
      #Quantiles Matching estimation
      mylocation <- as.numeric(median(obs))
      myscale <- as.numeric(IQR(obs, type=3)) / 2
      if(control$trace > 1)
      {
        cat("**\t default init\n")
        print(c(scale=as.numeric(myscale), location=as.numeric(mylocation)))
      }
    }
    
    denom <- myscale^2 + (obs - mylocation)^2
    scorestep <- c(mean(2*(obs - mylocation) / denom),
                   mean(1/myscale - 2*myscale / denom))
    LCE <- c(mylocation, myscale) + 2*myscale^2 * scorestep
    
    estimate <- c(location = LCE[1], scale = LCE[2])
    order  <- 0
    loglik <- NULL
    nbstep <- 1
  }
  
  if (distname == "weibull") 
  {
    if (initb)
    {
      shape <- init$shape
      scale <- init$scale
      #dweibull in stats pkg does not allow for a rate parameter
      rate <- 1/scale 
      if(control$trace > 1)
      {
        cat("**\t user supplied init\n")
        print(c(shape, scale))
      }
    }else
    {
      #if(is.null(delta))
      #  stop("wrong delta parameter")
      #ndelta <- as.integer(n^(delta))
      #idxsubsample <- sample(1:n, ndelta, replace=FALSE)
      #esttmp <- mledist(obs[idxsubsample], distr="weibull",...)$estimate
      #shape <- esttmp[1]
      #scale <- esttmp[2]
      #dweibull in stats pkg does not allow for a rate parameter
      #rate  <- 1/scale
      
      Y <- log(sort(obs))
      X <- log(-log(1-(1:n)/(n+1))) #voir aussi ASTM Procedure => utiliser ppoints()
      res <- as.numeric( coef( lm(Y~X) ) )
      shape <- 1/res[2]
      rate <- 1/exp(res[1])
      if(control$trace > 1)
      {
        cat("**\t default init\n")
        print(c(shape, scale))
      }
    }
    Score <- matrix(0, 2, 1)
    Score[1,1] <- mean(1/shape+ log(rate)+log(obs)-log(rate*obs)*(rate*obs)^shape)
    Score[2,1] <- mean(shape/rate-shape*(obs^shape)*rate^(shape-1))
    
    IFisher <- matrix(0,2,2)
    IFisher[1,1] <-  1/shape^2*(trigamma(1)+digamma(2)^2)
    IFisher[2,1] <-  1/rate*digamma(2) 
    IFisher[1,2] <-  1/rate*digamma(2)  
    IFisher[2,2] <-  shape^2/rate^2
    
    LCE <- c(shape,rate) + solve(IFisher, Score) 
    #dweibull in stats pkg does not allow for a rate parameter
    estimate <- c(shape = LCE[1], scale = 1/LCE[2])
    order <- 0
    loglik <- NULL
    nbstep <- 1
  } 
  
  if (distname == "t")
  {
    if (initb)
    {
      df<-init$df
      if(control$trace > 1)
      {
        cat("**\t user supplied init\n")
        print(df)
      }
    }else
    {
      if (delta <= 1 && delta > 1/2)
      {
        #subMLE initial guess estimator
        ndelta <- as.integer(n^(delta))
        idxsubsample <- 1:ndelta
        obstmp <- obs[idxsubsample]
        esttmp <- mledist(obstmp, "t", start=list("df"=10), control=con, ...)$estimate
        df <- esttmp["df"]
        
      }else if( delta < 1/2 && delta > 1/4 )
      {
        # Hill/MLE type initial guess estimator
        Y <- sort(obs, decreasing=TRUE)
        k <- floor(n^delta)
        Hilltmp <- mean(log(Y[1:k]))-log(Y[k+1])
        df <- 1/Hilltmp
      }else
      {
        stop("wrong value of delta")
      }
      if(control$trace > 1)
      {
        cat("**\t default init\n")
        print(df)
      }
    }
    
    #One Step 
    y <- obs 
    ysq <- y^2 #squared value
    Score <- 0.5*(digamma((df + 1)/2) - digamma(df / 2)) - 1/(2*df) + mean(- (1/2)*log(1+ysq/df) +(ysq*(df+1))/(2*df*(df + ysq)))
    
    IFisher <- (1/4)*(-trigamma((df+1)/2) + trigamma(df/2)) - (df+5)/(2*df*(df+1)*(df+3))
    
    LCE <- df + 1/IFisher * Score
    
    #Two Step 
    df <- LCE
    Score <- 0.5*(digamma((df + 1)/2) - digamma(df / 2)) - 1/(2*df) + mean(- (1/2)*log(1+ysq/df) +(ysq*(df+1))/(2*df*(df + ysq)))
    
    IFisher <- (1/4)*(-trigamma((df+1)/2) + trigamma(df/2)) - (df+5)/(2*df*(df+1)*(df+3))
    
    LCE <- df + 1/IFisher * Score
    estimate <- c(df = as.numeric(LCE)) 
    order <- 0
    loglik <- NULL
    nbstep <- 2
  }
  
  if (distname == "lst")
  {
    if (initb)
    {
      mu <- init$mu
      sigma <- init$sigma
      df <- init$df
      if(control$trace > 1)
      {
        cat("**\t user supplied init\n")
        print(c(df, sigma, mu))
      }
    }else
    {
      if (delta <= 1 && delta > 1/2)
      {
        #subMLE initial guess estimator
        ndelta <- as.integer(n^(delta))
        idxsubsample <- 1:ndelta
        obstmp <- obs[idxsubsample]
        startarg <- list("df"=10, "mu"=median(obstmp), "sigma"=IQR(obstmp)/2)
        esttmp <- mledist(obstmp, "lst", start=startarg, control=con, ...)$estimate
        
        df <- esttmp["df"]
        mu <- esttmp["mu"]
        sigma <- esttmp["sigma"]
        
      }else if( delta < 1/2 && delta > 1/4 )
      {
        # Hill/MLE type initial guess estimator
        mu <- median(obs)
        
        Y <- sort(obs-mu,decreasing=TRUE)
        k <- floor(n^delta)
        Hilltmp <- mean(log(Y[1:k]))-log(Y[k+1])
        df <- 1/Hilltmp
        
        minuslogvrais <- function(sigma)
        {
          ltmp <- dt(x=(obs-mu)/sigma, df=df, log=TRUE)
          lfinal <- ltmp[is.finite(ltmp)]
          res <- -sum(lfinal)+n*log(sigma)
          return(res)
        }
        
        #ressigma <- optim(IQR(obs)/2, minuslogvrais, control=con)
        ressigma <- try(optimize(minuslogvrais, lower=0, maximum=FALSE, tol=con$reltol,
                                 upper=quantile(obs, probs=3/4)))
        if(inherits(ressigma, "try-error"))
          sigma <- 1 #canonical value
        else
          sigma <- ressigma$minimum
        
      }else
      {
        stop("wrong value of delta")
      }
      if(control$trace > 1)
      {
        cat("**\t default init\n")
        print(c(df, sigma, mu))
      }
    }
    
    lst.score <- function(obs, mu, sigma, df)
    {
      y <- (obs-mu)/sigma #location-scaled obs
      ysq <- y^2 #squared value
      idx <- 1+ysq/df > 0
      Score <- matrix(0, 3, 1)
      Score[1,1] <- mean( (df+1)*y/(sigma*(df+ysq)) )
      Score[2,1] <- mean( -1/sigma + (df+1)*ysq/(sigma*(df+ysq)) )
      Score[3,1] <- 0.5*digamma((df + 1)/2) - 0.5*digamma(df / 2) - 1/(2*df) 
      Score[3,1] <- Score[3,1] + mean(- 0.5*log(1+ysq[idx]/df) + (ysq[idx]*(df+1))/(2*df*(df + ysq[idx])))
      Score
    }
    lst.Fisher <- function(sigma, df)
    {
      IFisher <- matrix(0, 3, 3)
      IFisher[1,1] <- (1+df)/((sigma^2)*(df+3))
      IFisher[2,2] <- 2*df/((sigma^2)*(df+3))
      IFisher[2,3] <- -2/(sigma*(df+1)*(df+3))
      IFisher[3,2] <- -2/(sigma*(df+1)*(df+3))
      IFisher[3,3] <- 0.25*(-trigamma((df+1)/2) + trigamma(df/2)) -(df+5)/(2*df*(df+1)*(df+3))
      IFisher
    }
    
    #One Step 
    Score <- lst.score(obs, mu, sigma, df)
    IFisher <- lst.Fisher(sigma, df)
    
    LCE1 <- c(mu, sigma, df) + solve(IFisher, Score) 
    
    #Two Step 
    mu <- LCE1[1]
    sigma <- LCE1[2]
    df <- LCE1[3]
    if(control$trace > 2)
    {
      cat("**\t after one-step\n")
      print(c(df, sigma, mu))
    }
    
    Score <- lst.score(obs, mu, sigma, df)
    IFisher <- lst.Fisher(sigma, df)
    
    LCE2 <- c(mu, sigma, df) + solve(IFisher, Score) 
    
    estimate <- c(mu = LCE2[1], sigma = LCE2[2], df= LCE2[3]) 
    order <- 0
    loglik <- NULL
    nbstep <- 2
  }
  
  if (distname == "pareto") 
  {
    if (initb)
    {
      shape <- init$shape
      scale <- init$scale
      if(control$trace > 1)
      {
        cat("**\t user supplied init\n")
        print(c(shape, scale))
      }
    }else
    {
      if(is.null(delta))
        stop("wrong delta parameter")
      ndelta <- as.integer(n^(delta))
      #idxsubsample <- sample(1:n, ndelta, replace=FALSE)
      idxsubsample <- 1:ndelta
      esttmp <- mledist(obs[idxsubsample], "pareto", control=con, ...)$estimate
      
      shape <- esttmp[1]
      scale <- esttmp[2]
      if(control$trace > 1)
      {
        cat("**\t default init\n")
        print(c(shape, scale))
      }
    }
    Score<- matrix(0,2,1)
    Score[1,1] <- mean(1/shape+log(scale)-log(scale+obs))
    Score[2,1] <- mean(shape/scale-(shape+1)*(obs+scale)^(-1))
    
    IFisher <- matrix(0,2,2)
    IFisher[1,1] <- 1/(shape^2)
    IFisher[2,1] <- -1/(scale*(shape+1))
    IFisher[1,2] <- -1/(scale*(shape+1)) 
    IFisher[2,2] <- shape/(scale^2*(shape+2))
    
    LCE <- c(shape, scale) + solve(IFisher, Score)
    estimate <- c(shape = LCE[1], scale = LCE[2])
    order <- 0
    loglik <- NULL
    nbstep <- 1
  }  
  list(estimate = estimate, convergence = 0, value = NULL, 
       hessian = NULL, optim.function = NULL, opt.meth = NULL, 
       fix.arg = NULL, fix.arg.fun = NULL, weights = NULL, 
       counts = NULL, optim.message = NULL, loglik = loglik, 
       method = "closed formula", order = order, memp = NULL,
       nbstep = nbstep)
}