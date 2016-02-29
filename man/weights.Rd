\name{weights}
\alias{weights.mvdareg}
\title{Extract Summary Information Pertaining to the Bootstrapped weights}
\description{Functions to extract \code{weights} bootstrap information from \code{mvdalab} objects.}
\usage{
\method{weights}{mvdareg}(object, ncomp = object$ncomp, conf = .95, ...)
}
\arguments{
\item{object}{ an mvdareg or mvdapaca object, i.e. \code{plsFit}. }
\item{ncomp}{ the number of components to include in the model (see below). }
\item{conf}{ for a bootstrapped model, the confidence level to use. }
\item{\dots}{ additional arguments.  Currently ignored. }
}
\details{
\code{weights} is used to extract a summary of the weights of a PLS.
If \code{ncomps} is missing (or is \code{NULL}), summaries for all regression estimates are returned.  Otherwise, if comps is given parameters for a model with only the requested component comps is returned. 

For mvdareg objects only, boostrap summaries provided are for actual regression weights, bootstrap percentiles, bootstrap mean, skewness, and bias.  These summaries can also be extracted using \code{weight.boots}
}
\value{
A weights object contains a data frame with columns:
\item{variable}{ variable names. }
\item{Actual}{ Actual loading estimate using all the data. }
\item{BCa percentiles}{ confidence intervals. }
\item{boot.mean}{ mean of the bootstrap. }
\item{skewness}{ skewness of the bootstrap distribution. }
\item{bias}{ estimate of bias w.r.t. the loading estimate. }
\item{Bootstrap Error}{ estimate of bootstrap standard error. }
\item{t value}{ approximate 't-value' based on the \code{Bootstrap Error}. }
\item{bias t value}{ approximate 'bias t-value' based on the \code{Bootstrap Error}. }
}
\author{Nelson Lee Afanador (\email{nelson.afanador@mvdalab.com})}
\seealso{\code{\link{weightsplot}}, \code{\link{weight.boots}}, \code{\link{weightsplot2D}}}
\references{
Davison, A.C. and Hinkley, D.V. (1997) Bootstrap Methods and Their Application. Cambridge University Press.

Efron, B. (1992) Jackknife-after-bootstrap standard errors and influence functions (with Discussion). Journal of the Royal Statistical Society, B, 54, 83:127.
}
\examples{
data(Penta)
mod1 <- plsFit(log.RAI ~., scale = TRUE, data = Penta[, -1], 
               ncomp = 3, contr = "contr.none", method = "bidiagpls", 
               validation = "oob")
weights(mod1, ncomp = 3, conf = .95)
}