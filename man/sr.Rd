\name{sr}
\alias{sr}
\alias{print.sr}
\alias{sr.error}
\alias{sr.modeled}
\title{Selectivity Ratio}
\description{This function calculates the Selectivity Ratio (\code{sr}) metric for an \code{mvdareg} object}
\usage{
sr(object, ncomps = object$ncomp)
}
\arguments{
\item{object}{ an mvdareg or mvdapaca object, i.e. \code{plsFit}. }
\item{ncomps}{ the number of components to include in the model (see below). }
}
\details{
\code{sr} is used to extract a summary of the significant multivariae correlation of a PLS model.

If \code{comps} is missing (or is \code{NULL}), summaries for all \code{sr} estimates are returned. Otherwise, if comps are given parameters for a model with only the requested component comps is returned. 
}
\value{
The output of \code{sr} is an \code{sr summary} detailing the following:
\item{sr}{ selectivity ratio statistic (\code{sr}).}
\item{p.value}{ p-value of the sr statistic. }
\item{f.value}{ f-value of the sr statistic. }
\item{Significant}{ Assessment of statistical significance. }
}
\author{Nelson Lee Afanador (\email{nelson.afanador@mvdalab.com})}
\references{
O.M. Kvalheim, T.V. Karstang, Interpretation of latent-variable regression models. Chemom. Intell. Lab. Syst., 7 (1989), pp. 39:51

O.M. Kvalheim, Interpretation of partial least squares regression models by means of target projection and selectivity ratio plots. J. Chemom., 24 (2010), pp. 496:504
}
\seealso{\code{\link{smc}}, \code{\link{vip}}}
\examples{
data(Penta)
mod1 <- plsFit(log.RAI ~., scale = TRUE, data = Penta[, -1], ncomp = 3, contr = "contr.none",
               method = "bidiagpls", validation = "oob")
sr(mod1)
plot(sr(mod1))
}