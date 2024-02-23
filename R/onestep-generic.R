onestep_generic <- function(obs, distname, ddistname, argddistname, control, init=init,...)
{
  n <- length(obs)
  
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
  
  if(!missing(init))
  {
    thetabar <- unlist(init)
    if(control$trace > 1)
    {
      cat("**\t user-supplied init\n")
      print(thetabar)
    }
  }else
  {
    #initial guess on a small subset of the original data
    delta <- control$delta
    if(is.null(delta))
      stop("wrong delta parameter")
    ndelta <- as.integer(n^(delta))
    idxsubsample <- sample(1:n, ndelta, replace=FALSE)
    
    #some optim control
    con <- list(trace = 0, fnscale = 1, maxit = 100L, abstol = -Inf, 
                reltol = sqrt(.Machine$double.eps), alpha = 1, beta = 0.5, 
                gamma = 2, REPORT = 10, warn.1d.NelderMead = TRUE, type = 1, 
                lmm = 5, factr = 1e+07, pgtol = 0, tmax = 10, temp = 10)
    con[names(control)] <- control
    con$trace <- max(con$trace-1, 0) #decrease trace
    
    resShortMLE <- try( mledist(data = obs[idxsubsample], distr = distname, 
                                control=con, ...) )
    if(inherits(resShortMLE, "try-error"))
    {
      return(list(estimate = NULL, convergence = 10, method = "numeric",
                  optim.message = "Initial guess cannot be estimated by MLE on a subsample"))
    }
    thetabar <- resShortMLE$estimate
    if(control$trace > 1)
    {
      cat("**\t sub-sample init\n")
      print(thetabar)
    }
  }
  
  #score definition
  Scorechap <- try( grad(fnobj, thetabar, obs=obs, ddistnam=ddistname) )
  if(inherits(Scorechap, "try-error"))
  {
    return(list(estimate = NULL, convergence = 10, method = "numeric",
                optim.message = "Gradient at initial guess cannot be estimated"))
  }
  
  #Hessian definition = opposite of Fisher information
  Ichap <- try( hessian(fnobj, thetabar, obs=obs, ddistnam=ddistname) )
  if(inherits(Ichap, "try-error"))
  {
    return(list(estimate = NULL, convergence = 10, method = "numeric",
                optim.message = "Hessian at initial guess cannot be estimated"))
  }
  
  #Compute a step of the Newton method
  step <- try( solve(Ichap, Scorechap) )
  if(inherits(step, "try-error"))
  {
    estimate <- NULL
    convergence <- 10
    optim.message <- "Newton step cannot be computed"
  }else
  {
    if(control$trace > 2)
    {
      cat("***\t step, Ichap, Scorechap\n")
      print(step)
      print(Ichap)
      print(Scorechap)
    }
    estimate <- thetabar - step
    convergence <- 0
    optim.message <- NULL
  }
  nbstep <- 1
  
  list(estimate = estimate, convergence = convergence, value = NULL, 
       hessian = NULL, optim.function = NULL, opt.meth = NULL, 
       fix.arg = NULL, fix.arg.fun = NULL, weights = NULL, 
       counts = NULL, optim.message = optim.message,
       method = "numeric", order = NULL, memp = NULL,
       nbstep = nbstep)
}