\name{coefsplot}
\alias{coefsplot}
\title{Graphical Summary Information Pertaining to the Regression Coefficients}
\description{Functions to extract regression coefficient bootstrap information from mvdalab objects.}
\usage{coefsplot(object, ncomp = object$ncomp, conf = 0.95)}
\arguments{
\item{object}{ an mvdareg object.  A fitted model. }
\item{ncomp}{ the number of components to include. }
\item{conf}{ for a bootstrapped model, the confidence level to use. }
}
\details{

\code{coefficients} is used to extract a graphical summary of the regression coefficients of a PLS model.

If \code{comps} is \code{missing} (or is \code{NULL}), a graphical summary for the nth component regression estimates are returned.  Otherwise, if \code{comps} is given parameters for a model with only the requested component \code{comps} is returned. 

Bootstrap graphcal summaries provided are when \code{method = oob}.

}
\author{Nelson Lee Afanador (\email{nelson.afanador@mvdalab.com})}
\examples{
data(Penta)
mod1 <- plsFit(log.RAI ~., scale = TRUE, data = Penta[, -1], 
        ncomp = 3, contr = "contr.none", method = "bidiagpls", 
        validation = "oob")
coefsplot(mod1, ncomp = 1:3)
}