benchonestep.replicate <- function(nsample, nbsimu, distr, methods=NULL,
                                     echo=FALSE, ncpus=1, ...) 
{
  rdistr <- paste0("r", distr)
  theoval <- list(...)
  if(is.null(methods))
    methods <- c("mme", "mle", "one")
    
  
  if(ncpus == 1)
  {
    bench1 <- function()
    {
      x <- do.call(rdistr, c(list(n=nsample), theoval))
      res <- benchonestep(data=x, distr=distr, methods=methods)
      parname <- rownames(res)[rownames(res) != "time"]
      error <- res[parname, ] - unlist(theoval)
      if(is.vector(error))
        error <- rbind(error)
      rownames(error) <- paste0("error-", parname)
      rbind(res, error)
    }
    res <- lapply(integer(nbsimu), 
                  function(n) bench1())
    res <- simplify2array(res)
    dimnames(res)[[3]] <- paste0("simu", 1:nbsimu)
    if(echo)
      cat("sample size", nsample, "\n")
  }else
  {
    bench1 <- 
      {
        ## force promises, so values get sent by parallel
        ## see line 154 of <boot>/R/bootfuns.q
        nsample; theoval; rdistr
        function()
        {
          x <- do.call(rdistr, c(list(n=nsample), theoval))
          res <- benchonestep(data=x, distr=distr, methods=methods)
          parname <- rownames(res)[rownames(res) != "time"]
          error <- res[parname, ] - unlist(theoval)
          if(is.vector(error))
            error <- rbind(error)
          rownames(error) <- paste0("error-", parname)
          rbind(res, error)
        }
      }
    cl <- parallel::makeCluster(ncpus)
    res <- parallel::parLapply(cl, integer(nbsimu), 
                               function(n) bench1())
    parallel::stopCluster(cl)
    res <- simplify2array(res)
    dimnames(res)[[3]] <- paste0("simu", 1:nbsimu)
  }
  res
}