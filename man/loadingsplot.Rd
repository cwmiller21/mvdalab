\name{loadingsplot}
\alias{loadingsplot}
\title{Graphical Summary Information Pertaining to the Loadings}
\description{Functions to extract graphical loadings information from \code{mvdareg} and \code{mvdapca} object.}
\usage{loadingsplot(object, ncomp = object$ncomp, conf = 0.95)}
\arguments{
\item{object}{ an \code{mvdareg} or \code{mvdapca} object.}
\item{ncomp}{ the number of components to include. }
\item{conf}{ for a bootstrapped model, the confidence level to use. }
}
\details{
\code{"loadingsplot"} is used to extract a graphical summary of the loadings of a PLS model.
If \code{"comps"} is missing (or is NULL), a graphical summary for the nth component estimates are returned.  Otherwise, if comps is given parameters for a model with only the requested component comps is returned. 

Bootstrap graphcal summaries provided are when \code{"method = oob"}
}
\author{Nelson Lee Afanador (\email{nelson.afanador@mvdalab.com})}
\seealso{\code{\link{loadings}}, \code{\link{loadings.boots}}, \code{\link{loadingsplot2D}}}
\examples{
data(Penta)
mod1 <- plsFit(log.RAI ~., scale = TRUE, data = Penta[, -1], 
               ncomp = 3, contr = "contr.none", method = "bidiagpls", 
               validation = "oob")
loadingsplot(mod1, ncomp = 1:3)
}