onestep <- function(data, distr, method, init, weights=NULL,
                    keepdata = TRUE, keepdata.nb=100, control=list(), ...)
{
  if (!is.character(distr))
  { 
    stop("distr must be a character string naming a distribution")
  }else
  {
    distname <- distr
  }	 
  
  distlist <- c("norm", "exp", "lnorm", "invgauss", "pois", "geom", "gamma","cauchy", 
                "nbinom", "beta", "unif", "logis", "weibull", "pareto", "chisq", "lst",
                "t")
  
  if(missing(method))
  {
    if(distname %in% distlist)
      method <- "closed formula"
    else 
      method <- "numeric"
  }
  method <- match.arg(method, c("closed formula", "numeric"))
  
  if (is.element(distname, c("binom", "nbinom", "geom", "hyper", "pois"))) 
    discrete <- TRUE
  else
    discrete <- FALSE
  
  #control parameters : param_t is the value at which the characteristic function is evaluated
  con <- list(delta = 0.9, 
              #optim control parameters
              trace=0, fnscale = 1, maxit = 100L, abstol = -Inf, 
              reltol = sqrt(.Machine$double.eps), alpha = 1, beta = 0.5, 
              gamma = 2, REPORT = 10, warn.1d.NelderMead = TRUE, type = 1, 
              lmm = 5, factr = 1e+07, pgtol = 0, tmax = 10, temp = 10)
  con[(namc <- names(control))] <- control
  
  if(length(con$delta) > 1)
    stop("delta control parameter is a scalar in [1/4, 1]")
  if(con$delta > 1 || con$delta < 1/4)
    stop("delta control parameter is a scalar in [1/4, 1]")
  
  if(is.element(distname, distlist) && method != "closed formula")
  {
    warning("argument 'method' could be set to 'closed formula'")
  }
  
  ddistname <- paste("d", distname, sep = "")
  if (!exists(ddistname, mode = "function")) 
    stop(paste("The ", ddistname, " function must be defined"))
  argddistname <- names(formals(ddistname))
  
  #check argument data
  if (!(is.vector(data) & is.numeric(data) & length(data)>1))
    stop("data must be a numeric vector of length greater than 1")
  obs <- data
  n <- length(obs)
  
  #check inconsistent initial parameters
  if(!missing(init))
  {
    hasnodefaultval <- sapply(formals(ddistname), is.name)
    arg_startfix <- checkparamlist(start.arg=init, fix.arg=NULL, 
      argdistname=argddistname, hasnodefaultval=hasnodefaultval)
  }
  
  #compute one-step estimator
  if (method == "closed formula") 
  {
    if(con$trace > 0)
      cat("*\t one step closed formula\n")
    res <- onestep_closedformula(obs=obs, distname=distname, control=con, init=init, ...)
    
  }else if(method == "numeric") 
  {  
    if(con$trace > 0)
      cat("*\t one step numeric\n")
    res <- onestep_generic(obs=obs, distname=distname, ddistname=ddistname, 
                           argddistname=argddistname, control=con, init=init,...)
    
  }else
    stop("wrong method")
  if(con$trace > 0)
  {
    cat("*\t one step result before processing output\n")
    print(res)
  }
  
  #compute components for a 'fitdist' object
  if(!is.null(res$hessian))
  {
    #check for NA values and invertible Hessian
    if(all(!is.na(res$hessian)) && qr(res$hessian)$rank == NCOL(res$hessian))
    {
      varcovar <- solve(res$hessian)
      sd <- sqrt(diag(varcovar))
      correl <- cov2cor(varcovar)
    }else
    {
      varcovar <- NA
      sd <- NA
      correl <- NA                            
    }
  }else
  {
    varcovar <- NA
    sd <- NA
    correl <- NA            
  }
  # Function to compute the log-likelihood to return
  if(is.null(res$loglik))
  {
    if(is.null(weights))
    {  
      if ("log" %in% argddistname)
      {
        loglik <- function(par, obs, ddistnam) 
          sum( do.call(ddistnam, c(list(obs), as.list(par), log=TRUE) ) )
      }else
      {
        loglik <- function(par, obs, ddistnam) 
          sum( log( do.call(ddistnam, c(list(obs), as.list(par)) ) ) )
      }
    }else #non NULL weight
    {
      if ("log" %in% argddistname)
      {
        loglik <- function(par, obs, ddistnam) 
          sum(weights * do.call(ddistnam, c(list(obs), as.list(par), log=TRUE) ) ) 
      }else
      {
        loglik <- function(par, obs, ddistnam) 
          sum(weights * log( do.call(ddistnam, c(list(obs), as.list(par)) ) ) )
      }
        
    }
    loglik <- loglik(res$estimate, obs, ddistname)
  }else
  {
    loglik <- res$loglik
  }
  
  npar <- length(res$estimate)
  aic <- -2*loglik+2*npar
  bic <- -2*loglik+log(n)*npar
  convergence <- res$convergence
  weights <- res$weights
  
  #encapsulate three dots arguments
  my3dots <- list(...)    
  if (length(my3dots) == 0) 
    my3dots <- NULL
  
  if(keepdata)
  {
    reslist <- list(estimate = res$estimate, method = method, sd = sd, cor = correl, 
                    vcov = varcovar, loglik = loglik, aic=aic, bic=bic, n = n, data=data,
                    distname = distname, fix.arg = NULL, fix.arg.fun = NULL, 
                    dots = my3dots, convergence = convergence, discrete = discrete, 
                    weights = res$Gts, nbstep = res$nbstep, delta=con$delta) 
  }else #just keep a sample set of all observations with min and max
  {
    n2keep <- min(keepdata.nb, n) - 2
    imin <- which.min(data)
    imax <- which.max(data)
    subdata <- data[sample((1:n)[-c(imin, imax)], size=n2keep, replace=FALSE)]
    subdata <- c(subdata, data[c(imin, imax)])
     
    reslist <- list(estimate = res$estimate, method = method, sd = sd, cor = correl, 
                  vcov = varcovar, loglik = loglik, aic=aic, bic=bic, n = n, data=subdata,
                  distname = distname, fix.arg = NULL, fix.arg.fun = NULL, 
                  dots = my3dots, convergence = convergence, discrete = discrete, 
                  weights = res$weights, nbstep = res$nbstep, delta=con$delta) 
  }
  if(con$trace > 1)
  {
    cat("**\t one step output\n")
    print(str(reslist))
  }
  
  class(reslist) <- c("onestep", "fitdist")
  
  return(reslist)
  
}
