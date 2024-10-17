\name{onestep}
\alias{onestep}
\title{
Executing Le Cam's one-step estimation
}
\description{
Executing Le Cam's one-step estimation based on Le Cam (1956) and Kamatani and Uchida (2015).
}
\usage{
onestep(data, distr, method, init, weights = NULL,  
  keepdata = TRUE, keepdata.nb=100, control=list(),  \dots)
}
\arguments{
  \item{data}{A numeric vector of length \code{n}}
  \item{distr}{A character string \code{"name"} naming a distribution for 
  which the corresponding density function dname and the 
  corresponding distribution function pname must be classically defined.}
  \item{method}{A character string coding for the fitting method:
    \code{"closed formula"} for explicit one-step and 
    \code{"numeric"} for numeric computation. The default method is the \code{"closed formula"}.}
  \item{init}{A named list for the user initial guess estimation.}
  \item{weights}{an optional vector of weights to compute the final likelihood. 
  Should be \code{NULL} or a numeric vector with strictly positive integers 
  (typically the number of occurences of each observation).}
  \item{keepdata}{a logical. If \code{TRUE},  dataset is returned, 
    otherwise only a sample subset is returned.}
  \item{keepdata.nb}{When \code{keepdata=FALSE},  the length (>1) of the subset returned.}
  \item{control}{a list of control parameters. Currently, 
  \code{param_t} is used when the characteristic function is needed,
  \code{delta} is used when the subsample of size \code{n^delta} is randomly selected
  for the initial guess in the generic Le Cam's one step method.}
  \item{\dots}{further arguments passe to \code{\link[fitdistrplus]{mledist}} in case it is used.}
}
\details{

The Le Cam one-step estimation procedure is based on an initial sequence of guess estimators and a Fisher scoring step or a single Newton step on the loglikelihood function. For the user,  the function onestep chooses automatically the best procedure to be used. The function OneStep presents internally several procedures depending on whether the sequence of initial guess estimators is in a closed form or not,  and on whether the score and the Fisher information matrix can be elicited in a closed form. "Closed formula" distributions are treated with explicit score and Fisher information matrix (or Hessian matrix). For all other distributions, if the density function is well defined, the numerical computation (\code{NumDeriv}) of the Newton step in Le Cam’s one-step is proposed with an initial sequence of guess estimators which is the sequence of maximum likelihood estimators computed on a subsample.
}
\value{
    \code{onestep} returns an object of class \code{"onestep"} inheriting from \code{"fitdist"}
    So,  it is a list with the following components:
    \item{estimate}{ the parameter estimates.}
    \item{method}{ the character string coding for the fitting method : 
        \code{"closed formula"} for closed-form MLE or closed-form one-step,  
        \code{"numeric"} for numeric computation of the one-step estimation.}
    \item{sd}{ the estimated standard errors,  \code{NA} if numerically not computable
        or \code{NULL} if not available.}
    \item{cor}{ the estimated correlation matrix,  \code{NA} if numerically not computable
        or \code{NULL} if not available.}
    \item{vcov}{ the estimated variance-covariance matrix,  \code{NULL} if not available.}
    \item{loglik}{ the log-likelihood.}
    \item{aic}{ the Akaike information criterion.}
    \item{bic}{ the the so-called BIC or SBC (Schwarz Bayesian criterion).}
    \item{n}{ the length of the data set.}
    \item{data}{ the data set.}
    \item{distname}{ the name of the distribution.}
    \item{dots}{ the list of  further arguments passed in \dots to be used .}
    \item{convergence}{ an integer code for the convergence:
        \code{0} indicates successful convergence (from explicit formula or not).
        \code{10} indicates an error.}    
    \item{discrete}{ the input argument or the automatic definition by the function to be passed
        to functions \code{\link[fitdistrplus]{gofstat}},  \code{\link[fitdistrplus]{plotdist}} 
        and \code{\link[fitdistrplus]{cdfcomp}}. }
    
    \item{weights}{the vector of weigths used in the estimation process or \code{NULL}.}    
    \item{nbstep}{the number of Newton step, 0 for closed-form MLE, 
      1 for one-step estimators and 2 for two-step estimators.}
    \item{delta}{delta parameter (used for sub-sample guess estimator).}  

    Generic functions inheriting from \code{"fitdist"} objects:
  \describe{       
    \item{\code{print}}{
    The print of a \code{"onestep"} object shows few traces about the fitting method and 
    the fitted distribution.
    }
    \item{\code{summary}}{
    The summary provides the parameter estimates of the fitted distribution,  the log-likelihood,  
    AIC and BIC statistics and when the maximum likelihood is used,  the standard errors of the 
    parameter estimates and the correlation matrix between parameter estimates.
    }
    \item{\code{plot}}{
    The plot of an object of class "onestep" returned by \code{fitdist} uses the function 
    \code{\link[fitdistrplus]{plotdist}}. An object of class "onestep" or a list of objects of class
    "onestep" corresponding to various fits using the same data set may also be plotted
    using a cdf plot (function \code{\link[fitdistrplus]{cdfcomp}}), 
    a density plot(function \code{\link[fitdistrplus]{denscomp}}), 
    a density Q-Q plot (function \code{\link[fitdistrplus]{qqcomp}}), 
    or a P-P plot (function \code{\link[fitdistrplus]{ppcomp}}).
    }
    \item{\code{logLik}}{ Extracts the estimated log-likelihood from the \code{"onestep"} object.
    }
    \item{\code{vcov}}{ Extracts the estimated var-covariance matrix from the \code{"onestep"} object.
    }
    \item{\code{coef}}{ Extracts the fitted coefficients from the \code{"onestep"} object.
    }
  }

}
\references{
L. Le Cam (1956). 
\emph{On the asymptotic theory of estimation and testing hypothesis}, 
In: Proceedings of 3rd Berkeley Symposium I, 355-368.

I.A. Koutrouvelis (1982). \emph{Estimation of Location and Scale in 
Cauchy Distributions Using the Empirical Characteristic Function}, 
Biometrika, 69(1), 205-213.

U. Grenander (1965). \emph{Some direct estimates of the mode}, 
Annals of Mathematical Statistics, 36, 131-138. 

K. Kamatani and M. Uchida (2015). \emph{Hybrid multi-step estimators 
for stochastic differential equations based on sampled data},
Stat Inference Stoch Process, 18(2), 177-204. 

Z.-S. Ye and N. Chen (2017). \emph{Closed-Form Estimators for the Gamma 
Distribution Derived From Likelihood Equations}, 
The American Statistician, 71(2), 177-181.


}
\author{
Alexandre Brouste, Christophe Dutang, Darel Noutsa Mieniedou 
}

\seealso{
See Also as \code{\link[fitdistrplus]{mledist}} and \code{\link[fitdistrplus]{fitdist}} in fitdistrplus.
}
\examples{
n <- 1000
set.seed(1234)

##1. Gamma

theta <- c(2, 3)
o.sample <- rgamma(n, shape=theta[1], rate=theta[2])

#Default method
onestep(o.sample, "gamma") 

#User initial sequence of guess estimator
# See : Ye and Chen (2017)

qtmp <- sum(o.sample*log(o.sample))-sum(log(o.sample))*mean(o.sample)
alphastar <- sum(o.sample)/qtmp
betastar <- qtmp/length(o.sample)
thetastar <- list(shape=alphastar,rate=1/betastar)  
  
onestep(o.sample, "gamma", init=thetastar) 

#Numerical method (for comparison)
onestep(o.sample, "gamma", method="numeric")

##2.Beta

theta <- c(0.5, 1.5)
o.sample <- rbeta(n, shape1=theta[1], shape2=theta[2])
onestep(o.sample, "beta")

##3. Cauchy

theta <- c(2, 3)
o.sample <- rcauchy(n, location=theta[1], scale=theta[2])
onestep(o.sample, "cauchy")

#User initial sequence of guess estimator
#See Koutrouvelis (1982).

t <- 1/4
Phi <- mean(exp(1i*t*o.sample))
S <- Re(Phi)
Z <- Im(Phi)
thetastar <- list(location=atan(Z/S)/t,scale=-log(sqrt(S^2+Z^2))/t)
onestep(o.sample, "cauchy", init=thetastar) 


##Chi2

theta <- 5
o.sample <- rchisq(n,df=theta)
onestep(o.sample,"chisq")

#User initial sequence of guess estimator
#See Grenander (1965).

p <- n^(2/7)
k <- floor(n^(3/5))
Dstar <- sort(o.sample)
Dk1 <- Dstar[(1+k):n]
Dk2 <- Dstar[1:(n-k)]
sigma.star <- 1/2*sum((Dk1+Dk2)*(Dk1-Dk2)^(-p))/sum((Dk1-Dk2)^(-p))+2

onestep(o.sample,"chisq",init=list(df=sigma.star))


#Negative Binomial

theta <- c(1, 5)
o.sample <- rnbinom(n, size=theta[1], mu=theta[2])
onestep(o.sample, "nbinom")


#Generic (dweibull2)

theta <- c(0.8, 3)
dweibull2 <- function(x, shape, scale, log=FALSE)
  dweibull(x = x, shape = shape, scale = scale, log = log)

o.sample <- rweibull(n, shape = theta[1], scale = 1/theta[2])
onestep(o.sample, "weibull2", method="numeric",  
  init=list(shape=1,  scale=1))
}
\keyword{ distribution }