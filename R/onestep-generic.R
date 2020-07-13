onestep_generic <- function(obs, distname, ddistname, argddistname, control, ...)
{
  n <- length(obs)
  delta <- control$delta
  if(is.null(delta))
    stop("wrong delta parameter")
  
  # generic definition of the log-likelihood => weights not taken into account
  if ("log" %in% argddistname) 
  {
    fnobj <- function(par, obs, ddistnam) 
    {
      -sum(do.call(ddistnam, c(list(obs), as.list(par), log = TRUE)))
    }
  }else 
  {
    fnobj <- function(par, obs, ddistnam) 
    {
      -sum(log(do.call(ddistnam, c(list(obs), as.list(par)))))
    }
  }
  
  #initial guess on a small subset of the original data
  ndelta <- as.integer(n^(delta))
  idxsubsample <- sample(1:n, ndelta, replace=FALSE)
  
  resShortMLE <- try( mledist(data = obs[idxsubsample], distr = distname, ...) )
  if(inherits(resShortMLE, "try-error"))
  {
    return(list(estimate = NULL, convergence = 10, method = "numeric",
                optim.message = "Initial guess cannot be estimated by MLE on a subsample"))
  }
  thetabar <- resShortMLE$estimate
  
  #score definition
  Scorechap <- try( grad(fnobj, thetabar, obs=obs, ddistnam=ddistname) )
  if(inherits(Scorechap, "try-error"))
  {
    return(list(estimate = NULL, convergence = 10, method = "numeric",
                optim.message = "Gradient at initial guess cannot be estimated"))
  }
  
  #Hessian definition
  Ichap <- try( hessian(fnobj, thetabar, obs=obs, ddistnam=ddistname) )
  if(inherits(Ichap, "try-error"))
  {
    return(list(estimate = NULL, convergence = 10, method = "numeric",
                optim.message = "Hessian at initial guess cannot be estimated"))
  }
  
  #onestep
  step <- try( solve(Ichap, Scorechap) )
  if(inherits(step, "try-error"))
  {
    estimate <- NULL
    convergence <- 10
    optim.message <- "Step cannot be computed"
  }else
  {
    estimate <- thetabar - step
    convergence <- 0
    optim.message <- NULL
  }
  
  list(estimate = estimate, convergence = convergence, value = NULL, 
            hessian = NULL, optim.function = NULL, opt.meth = NULL, 
            fix.arg = NULL, fix.arg.fun = NULL, weights = NULL, 
            counts = NULL, optim.message = optim.message,
            method = "numeric", order = NULL, memp = NULL)
}