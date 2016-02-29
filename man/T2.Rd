\name{T2}
\alias{T2}
\title{Generates a Hotelling's T2 Graph}
\description{Generates a Hotelling's T2 Graph both \code{mvdareg} and \code{mvdapca} objects.}
\usage{T2(object, ncomp = object$ncomp, phase = 1, conf = c(.95, .99))}
\arguments{
\item{object}{ an object of class \code{mvdareg} or \code{mvdapca}. }
\item{ncomp}{ the number of components to include in the calculation of Hotelling's T2. }
\item{phase}{ designates whether the confidence limits should reflect the current data frame, \code{phase = 1} or future observations, \code{phase = 2}. }
\item{conf}{ for a bootstrapped model, the confidence level to use. }
}
\details{
\code{T2} is used to generates a Hotelling's T2 graph both PLS and PCA models.
}
\value{
The output of \code{T2} is a graph of Hotelling's T2 and a data frame listing the T2 values.
}
\author{Nelson Lee Afanador (\email{nelson.afanador@mvdalab.com})}
\references{
Hotelling, H. (1931). "The generalization of Student's ratio". Annals of Mathematical Statistics 2 (3): 360:378.
}


\examples{
data(Penta)
mod1 <- plsFit(log.RAI ~., scale = TRUE, data = Penta[, -1], ncomp = 3, contr = "contr.none",
               method = "bidiagpls", validation = "oob")
T2(mod1, ncomp = 3)
}