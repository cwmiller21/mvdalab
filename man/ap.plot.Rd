\name{ap.plot}
\alias{ap.plot}
\title{Actual versus Predicted Plot and Residuals versus Predicted}
\description{This function provides the actual versus predicted and actual versus residuals plot as part of a model assessment}
\usage{ap.plot(object, ncomp = object$ncomp)}
\arguments{
\item{object}{ an object of class \code{mvdareg}, i.e., \code{plsFit}. }
\item{ncomp}{ number of components used in the model assessment }
}
\details{
This function provides the actual versus predicted and residuals versus predicted plot as part of model a assessment across the desired number of latent variables.
}
\value{
The output of \code{ap.plot} is a two facet graph for actual versus predicted and residuals versus predicted plots.
}
\seealso{\code{\link[mvdalab]{plsFit}}}
\author{Nelson Lee Afanador (\email{nelson.afanador@mvdalab.com})}
\examples{
mod1 <- plsFit(log.RAI ~., scale = TRUE, data = Penta[, -1], 
        ncomp = 3, contr = "contr.none", method = "bidiagpls", 
        validation = "oob")
ap.plot(mod1, ncomp = 3)
}