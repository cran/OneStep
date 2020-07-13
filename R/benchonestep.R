benchonestep <- function(data, distr, methods, init, weights=NULL,
                         ...) 
{
  restime <- numeric(length(methods))
  resestimate <- NULL
    
  for(i in 1:length(methods))
  {
    if(methods[i] == "mle")
    {
      restime[i] <- system.time(mymle <- mledist(data, distr)$estimate)[3]
      resestimate <- c(resestimate, list(mymle))
    }
    if(methods[i] == "mme")
    {
      restime[i] <- system.time(mymme <- mmedist(data, distr)$estimate)[3]
      resestimate <- c(resestimate, list(mymme))
    }
    if(methods[i] == "onestep")
    {
      restime[i] <- system.time(myos <- onestep(data, distr, method="closed form", init, weights=weights,
                                                 keepdata = FALSE, keepdata.nb=100, ...)$estimate)[3]
      resestimate <- c(resestimate, list(myos))
    }
  }
  resestimate <- simplify2array(resestimate)
  colnames(resestimate) <- methods
  
  
  return(rbind(resestimate, "time"=restime))
}
