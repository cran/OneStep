\name{benchonestep}
\alias{benchonestep}
\title{
Performing benchmark of one-step MLE against other methods
}
\description{
Performing benchmark of one-step MLE against other methods
}
\usage{
benchonestep(data, distr, methods, init, weights=NULL,\dots) 

}
\arguments{
  \item{data}{A numeric vector of length \code{n}}
  \item{distr}{	
A character string \code{"name"} naming a distribution for 
which the corresponding density function dname and the 
corresponding distribution function pname must be classically defined.}
  \item{methods}{A list of methods.}
  \item{init}{A character vector for the intial guess method.}
  \item{weights}{UNUSED : an optional vector of weights to be used in the fitting process. 
  Should be \code{NULL} or a numeric vector with strictly positive integers 
  (typically the number of occurences of each observation). 
  If non-\code{NULL}, weighted MLE is used, otherwise ordinary MLE.}
  \item{\dots}{unused}
}
\value{
    An array.

}
\references{
L. LeCam (1956). 
\emph{On the asymptotic theory of estimation and testing hypothesis}. 
In: Proceedings of 3rd Berkeley Symposium I, pages 355-368.
}
\author{
Alexandre Brouste, Darel Noutsa Mieniedou, Christophe Dutang
}

\keyword{ distribution }