\name{model.matrix}
\alias{model.matrix.mvdareg}
\title{\code{model.matrix} creates a design (or model) matrix.}
\description{This function returns the \code{model.matrix} of an \code{mvdareg} object.
}
\usage{
\method{model.matrix}{mvdareg}(object, ...)
}
\arguments{
\item{object}{ an \code{mvdareg} object }
\item{\dots}{ additional arguments.  Currently ignored. }
}
\details{
\code{"model.matrix.mvdareg"} is used to returns the \code{model.matrix} of an \code{mvdareg} object.
}
\value{
The design matrix for a \code{PLS} model with the specified formula and data.
}
\author{Nelson Lee Afanador (\email{nelson.afanador@mvdalab.com})}
\examples{
#PLS Model
mod1 <- plsFit(log.RAI ~., scale = TRUE, data = Penta[, -1], ncomp = 3, 
    contr = "contr.none", method = "bidiagpls", validation = "oob")
model.matrix(mod1)

}