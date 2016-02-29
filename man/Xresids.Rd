\name{Xresids}
\alias{Xresids}
\title{Generates a Graph of the X-residuals}
\description{Generates a graph of the X-residuals for both \code{mvdareg} and \code{mvdapca} objects.}
\usage{Xresids(object, ncomp = object$ncomp, conf = c(.95, .99), normalized = TRUE)}
\arguments{
\item{object}{ an object of class \code{mvdareg} or \code{mvdapca}. }
\item{ncomp}{ the number of components to include in the calculation of the X-residuals. }
\item{conf}{ for a bootstrapped model, the confidence level to use. }
\item{normalized}{ should residuals be normalized }
}
\details{
\code{Xresids} is used to generates a graph of the X-residuals for both PLS and PCA models.
}
\value{
The output of \code{Xresids} is a graph of X-residuals and a data frame listing the X-residuals values.
}
\author{Nelson Lee Afanador (\email{nelson.afanador@mvdalab.com})}
\references{
MacGregor, Process Monitoring and Diagnosis by Multiblock PLS Methods, May 1994 Vol. 40, No. 5 AIChE Journal.
}
\examples{
data(Penta)
mod1 <- plsFit(log.RAI ~., scale = TRUE, data = Penta[, -1], ncomp = 3, contr = "contr.none",
               method = "bidiagpls", validation = "oob")
Xresids(mod1, ncomp = 3)
}