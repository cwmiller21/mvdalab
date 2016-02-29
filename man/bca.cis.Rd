\name{bca.cis}
\alias{bca.cis}
\title{Bias-corrected and Accelerated Confidence Intervals}
\description{Computes bootstrap BCa confidence intervals for chosen parameters for PLS models fitted with \code{validation = "oob"}.}
\usage{bca.cis(object, conf = .95, type = c("coefficients", 
        "loadings", "weights"))}
\arguments{
\item{object}{ an object of class \code{"mvdareg"}, i.e. \code{plsFit}. }
\item{conf}{ desired confidence level }
\item{type}{ input parameter vector  }
}
\details{The function computes the bootstrap BCa confidence intervals for any fitted \code{mvdareg} model. 
Should be used in instances in which there is reason to suspectd the percentile intervals.  Results provided across all latent variables (LVs).  As such, it may be slow for models with a large number of LVs.
}
\value{
A bca.cis object contains component results for the following:
\item{ncomp}{ number of components in the model }
\item{variables}{ variable names }
\item{boot.mean}{ mean of the bootstrap }
\item{BCa percentiles}{ confidence intervals }
\item{proportional bias}{ calculated bias }
\item{skewness}{ skewness of the bootstrap distribution }
\item{a}{ acceleration contstant }
}
\author{Nelson Lee Afanador (\email{nelson.afanador@mvdalab.com})}
\references{
There are many references explaining the bootstrap and its implementation for confidence interval estimation. Among them are:

Davison, A.C. and Hinkley, D.V. (1997) Bootstrap Methods and Their Application. Cambridge University Press.

Efron, B. and Tibshirani, R. (1993) An Introduction to the Bootstrap. Chapman & Hall.

Hinkley, D.V. (1988) Bootstrap methods (with Discussion). Journal of the Royal Statistical Society, B, 50, 312:337, 355:370.
}

\seealso{\code{\link{plsFit}}, \code{\link{mvdaboot}}, \code{\link{boot.plots}}}
\examples{
data(Penta)
mod1 <- plsFit(log.RAI ~., scale = TRUE, data = Penta[, -1], 
        ncomp = 3, contr = "contr.none", method = "bidiagpls", 
        validation = "oob")
bca.cis(mod1, conf = .95, type = "coefficients")
bca.cis(mod1, conf = .95, type = "loadings")
bca.cis(mod1, conf = .95, type = "weights")
}