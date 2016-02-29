\name{mvrnorm.svd}
\alias{mvrnorm.svd}
\alias{mvrnormBase.svd}
\title{Simulate from a Multivariate Normal, Poisson, Exponential, or Skewed Distribution}
\description{Produces one or more samples from the specified multivariate  distribution.}
\usage{
mvrnorm.svd(n = 1, mu = NULL, Sigma = NULL, tol = 1e-06, empirical = FALSE, 
            Dist = "normal", skew = 5, skew.mean = 0, skew.sd = 1, 
            poisson.mean = 5)
}
\arguments{
\item{n}{ the number of samples required. }
\item{mu}{ a vector giving the means of the variables. }
\item{Sigma}{ a positive-definite symmetric matrix specifying the covariance matrix of the variables. }
\item{tol}{ tolerance (relative to largest variance) for numerical lack of positive-definiteness in Sigma. }
\item{empirical}{ logical. If true, \code{mu} and \code{Sigma} specify the empirical not population mean and covariance matrix. }
\item{Dist}{ desired distribution. }
\item{skew}{ amount of skew for skewed distributions. }
\item{skew.mean}{ mean for skewed distribution. }
\item{skew.sd}{ standard deviation for skewed distribution. }
\item{poisson.mean}{ mean for poisson distribution. }
}
\details{
\code{"mvrnorm.svd"} The matrix decomposition is done via \code{svd}
}
\author{Nelson Lee Afanador (\email{nelson.afanador@mvdalab.com})}
\examples{
Sigma <- matrix(c(1, .5, .5, .5, 1, .5, .5, .5, 1), 3, 3)
Means <- rep(0, 3)

Sim.dat.norm <- mvrnorm.svd(n = 1000, Means, Sigma, Dist = "normal")
plot(as.data.frame(Sim.dat.norm))

Sim.dat.pois <- mvrnorm.svd(n = 1000, Means, Sigma, Dist = "poisson")
plot(as.data.frame(Sim.dat.pois))

Sim.dat.exp <- mvrnorm.svd(n = 1000, Means, Sigma, Dist = "exp")
plot(as.data.frame(Sim.dat.exp))

Sim.dat.skew <- mvrnorm.svd(n = 1000, Means, Sigma, Dist = "skewnorm")
plot(as.data.frame(Sim.dat.skew))

}