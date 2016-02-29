\name{coef.mvdareg}
\alias{coef.mvdareg}
\title{Extract Information From a plsFit Model}
\description{Functions to extract information from \code{mvdalab} objects.}
\usage{
\method{coef}{mvdareg}(object, ncomp = object$ncomp, type = c("coefficients",
    "loadings", "weights", "y.loadings"), conf = .95, ...)
}
\arguments{
\item{object}{ an mvdareg object, i.e. a \code{plsFit}. }
\item{ncomp}{ the number of components to include in the model (see below). }
\item{type}{ specify model parameters to return. }
\item{conf}{ for a bootstrapped model, the confidence level to use. }
\item{\dots}{ additional arguments.  Currently ignored. }
}
\details{These are usually called through their generic functions \code{coef} and \code{residuals}, respectively.
\code{coef.mvdareg} is used to extract the regression coefficients, loadings, or weights of a PLS model.

If \code{comps} is missing (or is \code{NULL}), all parameter estimates are returned.
}
\value{
\item{coefficients}{ a named vector, or matrix, of coefficients. }
\item{loadings}{ a named vector, or matrix, of loadings. }
\item{weights}{ a named vector, or matrix, of weights. }
\item{y.loadings}{ a named vector, or matrix, of y.loadings. }
}
\author{Nelson Lee Afanador (\email{nelson.afanador@mvdalab.com})}
\seealso{\code{\link{coef}}, \code{\link{coefficients.boots}}, \code{\link{coefficients}}, 
\code{\link{loadings}}, \code{\link{loadings.boots}}, \code{\link{weights}}, 
\code{\link{weight.boots}}}
\examples{
mod1 <- plsFit(log.RAI ~., scale = TRUE, data = Penta[, -1], 
        ncomp = 3, contr = "contr.none", method = "bidiagpls", 
        validation = "oob")
coef(mod1, type = "coefficients", conf = .95)
}