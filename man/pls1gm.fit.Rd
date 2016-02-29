\name{pls1gm.fit}
\alias{pls1gm.fit}
\title{Modified Bidiag2 PLS}
\description{Fits a PLS model with the modified bidiagonalization algorithm for PLS1}
\usage{
pls1gm.fit(X, Y, ncomp, ...)
}
\arguments{
  \item{X}{ a matrix of observations. NAs and Infs are not allowed. }
  \item{Y}{ a vector. NAs and Infs are not allowed. }
  \item{ncomp}{ the number of components to include in the model (see below). }
  \item{\dots}{ additional arguments.  Currently ignored. }
}
\details{This function should not be called directly, but through the generic functions \code{plsFit} with the argument \code{method = "pls1gm"}. It implements the Bidiag2 scores algorithm.  
}
\value{
An object of class \code{mvdareg} is returned. The object contains all components returned by the underlying fit function. In addition, it contains the following components:
    \item{loadings}{ X loadings }
  \item{weights}{ weights }
  \item{D2.values}{ bidiag2 matrix }
  \item{iD2}{ inverse of bidiag2 matrix }
  \item{Ymean}{ mean of reponse variable }
  \item{Xmeans}{ mean of predictor variables }
  \item{coefficients}{ fitted values }
  \item{y.loadings}{ y-loadings }
  \item{scores}{ X scores }
  \item{R}{ orthogonal weights }
  \item{Y.values}{ scaled response values }
  \item{Yactual}{ actual response values }
  \item{fitted}{ fitted values }
  \item{residuals}{ residuals }
  \item{Xdata}{ X matrix }
  \item{iPreds}{ predicted values }
  \item{y.loadings2}{ scaled y-loadings }
  \item{ncomp}{ number of latent variables }
  \item{method}{ PLS algorithm used }
  \item{scale}{ scaling used }
  \item{validation}{ validation method }
  \item{call}{ model call }
  \item{terms}{ model terms }
  \item{model}{ fitted model }
}
\author{Nelson Lee Afanador (\email{nelson.afanador@mvdalab.com}), Thanh Tran (\email{thanh.tran@mvdalab.com})}
\references{
Indahl, Ulf G., (2014) The geometry of PLS1 explained properly: 10 key notes on mathematical properties of and some alternative algorithmic approaches to PLS1 modeling. Journal of Chemometrics, 28, 168:180.
}
\seealso{\code{\link{plsFit}}}