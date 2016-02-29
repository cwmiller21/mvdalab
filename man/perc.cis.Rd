\name{perc.cis}
\alias{perc.cis}
\title{Percentile Bootstrap Confidence Intervals}
\description{Computes percentile bootstrap confidence intervals for chosen parameters for \code{plsFit} models fitted with \code{validation = "oob"}}
\usage{
perc.cis(object, ncomp = object$ncomp, conf = 0.95, 
        type = c("coefficients", "loadings", "weights"))
}
\arguments{
\item{object}{ an object of class \code{"mvdareg"}, i.e., \code{plsFit} }
\item{ncomp}{ number of components to extract percentile intervals. }
\item{conf}{ confidence level. }
\item{type}{ input parameter vector. }
}
\details{The function fits computes the bootstrap percentile confidence intervals for any fitted \code{mvdareg} model. 

}
\value{
A perc.cis object contains component results for the following:
\item{ncomp}{ number of components in the model }
\item{variables}{ variable names }
\item{boot.mean}{ mean of the bootstrap }
\item{percentiles}{ confidence intervals }
}
\author{Nelson Lee Afanador (\email{nelson.afanador@mvdalab.com})}
\references{
There are many references explaining the bootstrap and its implementation for confidence interval estimation. Among them are:

Davison, A.C. and Hinkley, D.V. (1997) Bootstrap Methods and Their Application. Cambridge University Press.

Efron, B. and Tibshirani, R. (1993) An Introduction to the Bootstrap. Chapman & Hall.

Hinkley, D.V. (1988) Bootstrap methods (with Discussion). Journal of the Royal Statistical Society, B, 50, 312:337, 355:370.
}
\examples{
mod1 <- plsFit(log.RAI ~., scale = TRUE, data = Penta[, -1], 
        ncomp = 3, contr = "contr.none", method = "bidiagpls", 
        validation = "oob")
perc.cis(mod1, ncomp = 1:3, conf = .95, type = "coefficients")
}