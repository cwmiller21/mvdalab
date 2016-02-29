\name{acfplot}
\alias{acfplot}
\title{Plot of Auto-correlation Funcion}
\description{This function computes the autocorrelation function estimates for a selected parameter.}
\usage{acfplot(object, parm = NULL)}
\arguments{
\item{object}{ an object of class \code{mvdareg}, i.e., \code{plsFit}. }
\item{parm}{ a chosen predictor variable; if \code{NULL} a random predictor variable is chosen}
}
\details{
This function computes the autocorrelation function estimates for a selected parameter, via \code{acf}, and generates a graph that allows the analyst to assess the need for an autocorrelation adjustment in the \code{\link{smc}}.
}
\author{Nelson Lee Afanador (\email{nelson.afanador@mvdalab.com})}
\references{
This function is built using the \code{acf} function in the \pkg{stats} R package.

Venables, W. N. and Ripley, B. D. (2002) Modern Applied Statistics with S. Fourth Edition. Springer-Verlag.
}
\seealso{\code{\link{smc}}, \code{\link{smc.acfTest}}}
\examples{
data(Penta)
mod1 <- plsFit(log.RAI ~., scale = TRUE, data = Penta[, -1], 
        ncomp = 3, contr = "contr.none", method = "bidiagpls", 
        validation = "oob")
acfplot(mod1, parm = NULL)
}