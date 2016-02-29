\name{weightsplot}
\alias{weightsplot}
\title{Extract Graphical Summary Information Pertaining to the Weights}
\description{Functions to extract regression coefficient bootstrap information from \code{mvdalab} objects.}
\usage{weightsplot(object, ncomp = object$ncomp, conf = .95)}
\arguments{
\item{object}{ an mvdareg object, i.e. \code{plsFit} }
\item{ncomp}{ the number of components to include. }
\item{conf}{ for a bootstrapped model, the confidence level to use. }
}
\details{

\code{weightsplot} is used to extract a graphical summary of the weights of a PLS model.

If \code{comps} is missing (or is \code{NULL}), a graphical summary for the nth component regression estimates are returned.  Otherwise, if comps is given parameters for a model with only the requested component \code{comps} is returned. 

Boostrap graphcal summaries provided are when \code{method = oob}.
}
\author{Nelson Lee Afanador (\email{nelson.afanador@mvdalab.com})}
\examples{
data(Penta)
mod1 <- plsFit(log.RAI ~., scale = TRUE, data = Penta[, -1], ncomp = 3, contr = "contr.none",
               method = "bidiagpls", validation = "oob")
weightsplot(mod1, ncomp = 1:3)
}