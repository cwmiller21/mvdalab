\name{y.loadings.boots}
\alias{y.loadings.boots}
\title{Extract Summary Information Pertaining to the y-loadings}
\description{Functions to extract the y-loadings from mvdareg and mvdapca objects.}
\usage{y.loadings.boots(object, ncomp = object$ncomp, conf = 0.95)}
\arguments{
\item{object}{ an \code{mvdareg} or \code{mvdapaca} object, i.e. \code{plsFit}. }
\item{ncomp}{ the number of components to include in the model (see below). }
\item{conf}{ for a bootstrapped model, the confidence level to use. }
}
\details{

\code{y.loadings.boots} is used to extract a summary of the y-loadings from a PLS or PCA model.

If \code{comps} is missing (or is \code{NULL}), summaries for all regression estimates are returned.  Otherwise, if \code{comps} is provided the requested component comps are returned. 

For \code{mvdareg} objects only, boostrap summaries provided are for actual regression \code{y.loadings}, bootstrap percentiles, bootstrap mean, skewness, and bias.  These summaries can also be extracted using \code{y.loadings.boots}

}
\author{Nelson Lee Afanador (\email{nelson.afanador@mvdalab.com})}
\examples{
data(Penta)
mod1 <- plsFit(log.RAI ~., scale = TRUE, data = Penta[, -1], ncomp = 3, contr = "contr.none",
               method = "bidiagpls", validation = "oob")
y.loadings(mod1)
y.loadings.boots(mod1)
}