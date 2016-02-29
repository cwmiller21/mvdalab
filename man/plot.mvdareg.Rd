\name{plot.mvdareg}
\alias{plot.mvdareg}
\title{General plotting function for \code{mvdareg} and \code{mvdapaca} objects.}
\description{A general plotting function for a \code{mvdareg} and \code{mvdapca} objects.}
\usage{
\method{plot}{mvdareg}(x, plottype = c("PE", "scoresplot", "loadingsplot", 
                    "loadingsplot2D", "T2", "Xresids", "coefsplot", "ap.plot", 
                    "weightsplot", "weightsplot2D", "acfplot"), ...)
}
\arguments{
\item{x}{ an object of class \code{"mvdareg"}, i.e., a fitted model. }
\item{plottype}{ the desired plot from an object of class \code{"mvdareg"} }
\item{\dots}{ additional arguments.  Currently ignored. }
}
\details{
The following plotting functions are supported:

\code{PE}, \code{scoreplot}, \code{loadingsplot}, \code{loadingsplot2D}, \code{T2}, \code{Xresids}, \code{coefsplot}, \code{ap.plot}, \code{weightsplot}, \code{weightsplot2D}, \code{acfplot}
}
\author{Nelson Lee Afanador (\email{nelson.afanador@mvdalab.com})}

\examples{
mod1 <- plsFit(log.RAI ~., scale = TRUE, data = Penta[, -1], 
        ncomp = 3, contr = "contr.none", method = "bidiagpls", 
        validation = "oob")
plot(mod1, plottype = "PE")
plot(mod1, plottype = "scoresplot")
plot(mod1, plottype = "loadingsplot", ncomp = 1:2)
plot(mod1, plottype = "loadingsplot2D")
plot(mod1, plottype = "T2", ncomp = 2, phase = 1, conf = c(.95, .99))
plot(mod1, plottype = "Xresids")
plot(mod1, plottype = "coefsplot")
plot(mod1, plottype = "ap.plot")
plot(mod1, plottype = "weightsplot")
plot(mod1, plottype = "weightsplot2D")
plot(mod1, plottype = "acfplot", parm = "L2")
}