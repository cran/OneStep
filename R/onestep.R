onestep <- function(data, distr, method, init, weights=NULL,
                    keepdata = TRUE, keepdata.nb=100, control=list(), ...)
{
  
  if (!is.character(distr))
  { 
    stop("distr must be a character string naming a distribution")
  }else{
    distname <- distr
  }	 
  
  distlist <- c("norm", "exp", "lnorm", "invgauss", "pois", "geom", "gamma","cauchy", 
                "nbinom", "beta", "unif", "logis","weibull","pareto")
  if(missing(method))
  {
    if(distname %in% distlist)
      method <- "closed formula"
    else 
      method <- "numeric"
  }
    
  if (is.element(distname, c("binom", "nbinom", "geom", "hyper", "pois"))) 
    discrete <- TRUE
  else
    discrete <- FALSE
  
  method <- match.arg(method, c("closed formula", "numeric"))
  #control parameters : param_t is the value at which the characteristic function is evaluated
  con <- list(param_t = 0.3, delta=0.8)
  con[(namc <- names(control))] <- control
  
  if (is.element(distname, distlist) && method != "closed formula")
  {
    warning("argument 'method' could be set to 'closed formula'")
  }
  
  ddistname <- paste("d", distname, sep = "")
  argddistname <- names(formals(ddistname))
  if (!exists(ddistname, mode = "function")) 
    stop(paste("The ", ddistname, " function must be defined"))
  
  obs <- data
  n <- length(obs)
  
  #compute one-step estimator
  if (method == "closed formula") 
  {
    res <- onestep_closedformula(data=obs, distname=distname, control=con, ...)
    
  }else if(method == "numeric") 
  {  
    res <- onestep_generic(obs=obs, distname=distname, ddistname=ddistname, 
                           argddistname=argddistname, control=con, ...)
    
  }else
    stop("wrong method")
  
  
  #compute components for a 'fitdist' object
  if(!is.null(res$hessian)){
    #check for NA values and invertible Hessian
    if(all(!is.na(res$hessian)) && qr(res$hessian)$rank == NCOL(res$hessian)){
      varcovar <- solve(res$hessian)
      sd <- sqrt(diag(varcovar))
      correl <- cov2cor(varcovar)
    }else{
      varcovar <- NA
      sd <- NA
      correl <- NA                            
    }
  }else{
    varcovar <- NA
    sd <- NA
    correl <- NA            
  }
  # Function to calculate the loglikelihood to return
  if(is.null(weights))
  {  
    loglik <- function(par, obs, ddistnam) 
      sum(log(do.call(ddistnam, c(list(obs), as.list(par)) ) ) )
  }else
  {
    loglik <- function(par, obs, ddistnam) 
      sum(weights * log(do.call(ddistnam, c(list(obs), as.list(par)) ) ) )
  }
  
  loglik <- loglik(res$estimate, obs, ddistname)
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
                    weights = res$Gts) 
  }else#just keep a sample set of all observations
  {
    n2keep <- min(keepdata.nb, n)-2
    imin <- which.min(data)
    imax <- which.max(data)
    subdata <- data[sample((1:n)[-c(imin, imax)], size=n2keep, replace=FALSE)]
    subdata <- c(subdata, data[c(imin, imax)])
     
    reslist <- list(estimate = res$estimate, method = method, sd = sd, cor = correl, 
                  vcov = varcovar, loglik = loglik, aic=aic, bic=bic, n = n, data=subdata,
                  distname = distname, fix.arg = NULL, fix.arg.fun = NULL, 
                  dots = my3dots, convergence = convergence, discrete = discrete, 
                  weights = res$weights) 
  }
  class(reslist) <- c("onestep", "fitdist")
  
  return(reslist)
  
}
