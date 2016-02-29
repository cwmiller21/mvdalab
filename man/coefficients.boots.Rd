\name{coefficients.boots}
\alias{coefficients.boots}
\title{BCa Summaries for the coefficient of an mvdareg object}
\description{Computes bootstrap BCa confidence intervals for regression coefficients, along with expanded bootstrap summaries.}
\usage{coefficients.boots(object, ncomp = object$ncomp, conf = 0.95)}
\arguments{
\item{object}{ an object of class \code{mvdareg}, i.e., a \code{plsFit}. }
\item{ncomp}{ number of components in the model  }
\item{conf}{ desired confidence level }
}
\details{The function computes the bootstrap BCa confidence intervals for fitted \code{mvdareg} models where \code{valiation = "oob"}. 
Should be used in instances in which there is reason to suspectd the percentile intervals.  Results provided across all latent variables or for specific latent variables via \code{ncomp}.  
}
\value{
A coefficients.boots object contains component results for the following:
\item{variable}{ variable names }
\item{actual}{ Actual loading estimate using all the data }
\item{BCa percentiles}{ confidence intervals }
\item{boot.mean}{ mean of the bootstrap }
\item{skewness}{ skewness of the bootstrap distribution }
\item{bias}{ estimate of bias w.r.t. the loading estimate }
\item{Bootstrap Error}{ estimate of bootstrap standard error }
\item{t value}{ approximate 't-value' based on the \code{Bootstrap Error} }
\item{bias t value}{ approximate 'bias t-value' based on the \code{Bootstrap Error} }
}
\author{Nelson Lee Afanador (\email{nelson.afanador@mvdalab.com})}
\seealso{\code{\link{coef}}, \code{\link{coefficients}}, 
\code{\link{coefsplot}}, \code{\link{coefficients}}}
\references{

There are many references explaining the bootstrap. Among them are:

Davison, A.C. and Hinkley, D.V. (1997) Bootstrap Methods and Their Application. Cambridge University Press.

Efron, B. and Tibshirani, R. (1993) An Introduction to the Bootstrap. Chapman & Hall.

Hinkley, D.V. (1988) Bootstrap methods (with Discussion). Journal of the Royal Statistical Society, B, 50, 312:337, 355:370.

}
\examples{
mod1 <- plsFit(log.RAI ~., scale = TRUE, data = Penta[, -1], 
        ncomp = 3, contr = "contr.none", method = "bidiagpls", 
        validation = "oob")
coefficients.boots(mod1, ncomp = 3, conf = .95)
}