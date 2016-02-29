\name{plot.R2s}
\alias{plot.R2s}
\title{Plot of R2}
\description{Plots for the cross-validated R2 (CVR2), explained variance in the predictor variables (R2X), and the reponse (R2Y).
}
\usage{
\method{plot}{R2s}(x, ...)
}
\arguments{
\item{x}{An \code{R2s} object }
\item{\dots}{ additional arguments.  Currently ignored. }
}
\details{
\code{plot.R2s} is used to generates the graph of the cross-validated R2 (CVR2), explained variance in the predictor variables (R2X), and the reponse (R2Y) for PLS models.

}
\value{
The output of \code{plot.R2s} is a graph of the stated explained variance summary.
}
\author{Thanh Tran (\email{thanh.tran@mvdalab.com})}
\examples{
data(Penta)
mod1 <- plsFit(log.RAI ~., scale = TRUE, data = Penta[, -1], ncomp = 3, contr = "contr.none",
               method = "bidiagpls", validation = "oob")
R2s. <- R2s(mod1)
plot(R2s.)

mod2 <- plsFit(log.RAI ~., scale = TRUE, data = Penta[, -1], ncomp = 3, contr = "contr.none",
               method = "bidiagpls", validation = "loo")
R2s.. <- R2s(mod2)
plot(R2s..)
}