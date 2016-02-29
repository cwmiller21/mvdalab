\name{plot.smc}
\alias{plot.smc}
\title{Plotting function for Significant Multivariate Correlation}
\description{This function generates a plot an object of class \code{smc}.
}
\usage{
\method{plot}{smc}(x, variables = "all", ...)
}
\arguments{
\item{x}{\code{smc} object. }
\item{variables}{ the number of variables to include the graph output. }
\item{\dots}{ additional arguments.  Currently ignored. }
}
\details{

\code{plot.smc} is used to generates the graph of the significant multivariate correlation from \code{smc} objects.
}
\value{
The output of \code{plot.smc} is a graph of the significant multivariate correlation for the specified observation(s).
}
\author{Nelson Lee Afanador (\email{nelson.afanador@mvdalab.com})}

\examples{
data(Penta)
mod1 <- plsFit(log.RAI ~., scale = TRUE, data = Penta[, -1], ncomp = 3, contr = "contr.none",
               method = "bidiagpls", validation = "oob")
smc(mod1)
plot(smc(mod1))
}