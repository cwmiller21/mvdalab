\name{print.mvdalab}
\alias{print.mvdareg}
\title{Print Methods for mvdalab Objects}
\description{Summary and print methods for mvdalab objects.}
\usage{
\method{print}{mvdareg}(x, ...)
}
\arguments{
\item{x}{ an mvdalab object }
\item{\dots}{ additional arguments.  Currently ignored. }
}
\details{
\code{print.mvdalab} Is a generic function used to print mvdalab objects, such as \code{print.empca} for \code{imputeEM}, \code{print.mvdapca} for \code{mvdapca} objects, and \code{summary.mvdareg} for \code{mvdareg} objects.
}
\author{Nelson Lee Afanador (\email{nelson.afanador@mvdalab.com})}
\examples{
data(Penta)
mod1 <- plsFit(log.RAI ~., scale = TRUE, data = Penta[, -1], ncomp = 3, contr = "contr.none",
               method = "bidiagpls", validation = "oob")
print(mod1, ncomp = 3)
summary(mod1, ncomp = 3)
}