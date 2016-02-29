\name{XresidualContrib}
\alias{XresidualContrib}
\title{Generates the squared prediction error contributions and contribution plot}
\description{Generates the squared prediction error (SPE) contributions and graph both \code{mvdareg} and \code{mvdapca} objects.}
\usage{XresidualContrib(object, ncomp = object$ncomp, obs1 = 1)}
\arguments{
\item{object}{ an object of class \code{mvdareg} or \code{mvdapca}. }
\item{ncomp}{ the number of components to include in the SPE calculation. }
\item{obs1}{ the observaion in SPE assessment. }

}
\details{
\code{XresidualContrib} is used to generates the squared prediction error (SPE) contributions and graph for both PLS and PCA models.  Only one observation at a time is supported.


}
\value{
The output of \code{XresidualContrib} is a matrix of score contributions for a specified observation and the corresponding graph.
}
\author{Nelson Lee Afanador (\email{nelson.afanador@mvdalab.com})}
\references{
MacGregor, Process Monitoring and Diagnosis by Multiblock PLS Methods, May 1994 Vol. 40, No. 5 AIChE Journal
}


\examples{
#PLS Model
#PLS Model
mod1 <- plsFit(log.RAI ~., scale = TRUE, data = Penta[, -1], ncomp = 3, contr = "contr.none",
               method = "bidiagpls", validation = "oob")
XresidualContrib(mod1, ncomp = 3, obs1 = 3)
XresidualContrib(mod1, ncomp = 3, obs1 = 5)


#PCA Model
pc1 <- pcaFit(Penta[, -1])
XresidualContrib(mod1, ncomp = 3, obs1 = 3)
XresidualContrib(mod1, ncomp = 3, obs1 = 5)
}