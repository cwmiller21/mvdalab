\name{R2s}
\alias{R2s}
\alias{print.R2s}
\title{Cross-validated R2, R2 for X, and R2 for Y for PLS models}
\description{Functions to report the cross-validated R2 (CVR2), explained variance in the predictor variables (R2X), and the reponse (R2Y) for PLS models.}
\usage{R2s(object)}
\arguments{
\item{object}{ an mvdareg object, i.e., \code{plsFit}. }
}
\details{
\code{R2s} is used to extract a summary of the cross-validated R2 (CVR2), explained variance in the predictor variables (R2X), and the reponse (R2Y) for PLS models.

}
\author{Nelson Lee Afanador (\email{nelson.afanador@mvdalab.com})}
\examples{
data(Penta)
mod1 <- plsFit(log.RAI ~., scale = TRUE, data = Penta[, -1], ncomp = 3, contr = "contr.none",
               method = "bidiagpls", validation = "oob")
R2s(mod1)
plot(R2s(mod1))
}