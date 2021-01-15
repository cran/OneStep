benchonestep <- function(data, distr, methods, init, weights=NULL,
                         ...) 
{
  
  #check argument data
  if (!(is.vector(data) & is.numeric(data) & length(data)>1))
    stop("data must be a numeric vector of length greater than 1")
  #check argument distr
  if (!is.character(distr))
  { 
    stop("distr must be a character string naming a distribution")
  }else{
    distname <- distr
  }	 
  #check argument methods
  methods <- sapply(methods, function(x) 
    match.arg(x, c("mme", "mle", "onestep")))
  
  restime <- numeric(length(methods))
  resestimate <- NULL
    
  for(i in 1:length(methods))
  {
    if(methods[i] == "mle")
    {
      timebeg <- proc.time()
      mymle <- mledist(data, distr, weights=weights)$estimate
      timeend <- proc.time()
      restime[i] <- timeend[3] - timebeg[3]
      resestimate <- c(resestimate, list(mymle))
    }
    if(methods[i] == "mme")
    {
      timebeg <- proc.time()
      mymme <- mmedist(data, distr, weights=weights)$estimate
      timeend <- proc.time()
      restime[i] <- timeend[3] - timebeg[3]
      resestimate <- c(resestimate, list(mymme))
    }
    if(methods[i] == "onestep")
    {
      timebeg <- proc.time()
      myos <- onestep(data, distr, method="closed form", init, 
                      weights=weights, keepdata = FALSE, keepdata.nb=100, ...)$estimate
      timeend <- proc.time()
      restime[i] <- timeend[3] - timebeg[3]
      resestimate <- c(resestimate, list(myos))
    }
  }
  resestimate <- simplify2array(resestimate)
  if(is.vector(resestimate))
  {
    resestimate <- rbind(resestimate)
    rownames(resestimate) <- colnames(resestimate)[1]
  }
  colnames(resestimate) <- methods
  
  return(rbind(resestimate, "time"=restime))
}
