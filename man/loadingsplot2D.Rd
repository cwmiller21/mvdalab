\name{loadingsplot2D}
\alias{loadingsplot2D}
\title{2-Dimensionsl Graphical Summary Information Pertaining to the Loadings of a PLS or PCA Analysis}
\description{Functions to extract 2D graphical loadings information from mvdalab objects.}
\usage{
loadingsplot2D(object, comps = c(1, 2))
}
\arguments{
\item{object}{ an \code{mvdareg} or \code{mvdapca} object. }
\item{comps}{ a vector or length 2 corresponding to the number of components to include. }
}
\details{
\code{loadingsplot2D} is used to extract a graphical summary of the loadings of a PLS model.
If \code{comp} is missing (or is NULL), a graphical summary for the 1st and 2nd componentsare returned.  
}
\author{Nelson Lee Afanador (\email{nelson.afanador@mvdalab.com})}
\seealso{\code{\link{coefficientsplot2D}}, \code{\link{weightsplot2D}}}
\examples{
data(Penta)
mod1 <- plsFit(log.RAI ~., scale = TRUE, data = Penta[, -1], ncomp = 3, contr = "contr.none",
               method = "bidiagpls", validation = "oob")
loadingsplot2D(mod1, comp = c(1, 2))

data(iris)
pc1 <- pcaFit(iris)
loadingsplot2D(pc1, comp = c(1, 2))
}