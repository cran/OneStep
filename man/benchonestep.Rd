\name{benchonestep}
\alias{benchonestep}
\alias{benchonestep.replicate}
\title{
Performing benchmark of one-step MLE against other methods
}
\description{
\code{benchonestep} performs a benchmark of one-step MLE against other methods on a given dataset.
\code{benchonestep.replicate} repeats several times the procedure: data random generation and benchmark through \code{benchonestep}.
}
\usage{
benchonestep(data, distr, methods, init, weights=NULL,\dots) 
benchonestep.replicate(nsample, nbsimu, distr, methods=NULL, echo=FALSE, ncpus=1, \dots) 
}
\arguments{
  \item{data}{A numeric vector of length \code{n}}
  \item{distr}{A character string \code{"name"} naming a distribution for 
  which the corresponding density function dname and the 
  corresponding distribution function pname must be classically defined.}
  \item{methods}{A vector of methods: character among
    \code{"mme"}, \code{"mle"}, \code{"onestep"} (can be abbreviated).}
  \item{init}{A named list for the intial guess method.}
  \item{weights}{An optional vector of weights to be used in the fitting process. 
  Should be \code{NULL} or a numeric vector with strictly positive integers 
  (typically the number of occurences of each observation). 
  If non-\code{NULL}, weighted MLE is used, otherwise ordinary MLE.}
  \item{\dots}{unused for \code{benchonestep}; true parameters passed to \code{rdistr}
  for \code{benchonestep.replicate}}
  \item{nsample}{a numeric for the sample size.}
  \item{nbsimu}{a numeric for the replication number.}
  \item{echo}{a logical to display or not some traces of benchmarking.}
  \item{ncpus}{Number of processes to be used in parallel operation: 
  typically one would fix it to the number of available CPUs.}
}
\value{
    A matrix with estimate and time in seconds per method for  \code{benchonestep};
    an array with estimates, times, errors in seconds per method, per simulation for  \code{benchonestep.replicate}.

}
\references{
L. LeCam (1956). 
\emph{On the asymptotic theory of estimation and testing hypothesis}. 
In: Proceedings of 3rd Berkeley Symposium I, pages 355-368.
}
\author{
Alexandre Brouste, Darel Noutsa Mieniedou, Christophe Dutang
}
\examples{
n <- 1000
set.seed(1234)

x <- rbeta(n, 3, 2)
benchonestep(x, "beta", c("mle", "one")) 
}
\keyword{ distribution }