\name{smc}
\alias{smc}
\alias{print.smc}
\alias{smc.error}
\alias{smc.modeled}
\title{Significant Multivariate Correlation}
\description{This function calculates the significant multivariate correlation (\code{smc}) metric for an \code{mvdareg} object}
\usage{
smc(object, ncomps = object$ncomp, corrected = F)
}
\arguments{
\item{object}{ an mvdareg or mvdapaca object, i.e. \code{plsFit}. }
\item{ncomps}{ the number of components to include in the model (see below). }
\item{corrected}{ whether there should be a correction of 1st order auto-correlation in the residuals. }
}
\details{
\code{smc} is used to extract a summary of the significant multivariae correlation of a PLS model.

If \code{comps} is missing (or is \code{NULL}), summaries for all \code{smc} estimates are returned. Otherwise, if comps are given parameters for a model with only the requested component comps is returned. 
}
\value{
The output of \code{smc} is an \code{smc summary} detailing the following:

\item{smc}{ significant multivariate correlation statistic (\code{smc}).}
\item{p.value}{ p-value of the smc statistic. }
\item{f.value}{ f-value of the smc statistic. }
\item{Significant}{ Assessment of statistical significance. }

}
\author{Nelson Lee Afanador (\email{nelson.afanador@mvdalab.com})}
\references{

Thanh N. Tran, Nelson Lee Afanador, Lutgarde M.C. Buydens, Lionel Blanchet, Interpretation of variable importance in Partial Least Squares with Significance Multivariate Correlation (sMC). Chemom. Intell. Lab. Syst. 2014; 138: 153:160.

Nelson Lee Afanador, Thanh N. Tran, Lionel Blanchet, Lutgarde M.C. Buydens, Variable importance in PLS in the presence of autocorrelated data - Case studies in manufacturing processes. Chemom. Intell. Lab. Syst. 2014; 139: 139:145.

}
\seealso{\code{\link{smc.acfTest}}, \code{\link{sr}}, \code{\link{vip}}}
\examples{
data(Penta)
mod1 <- plsFit(log.RAI ~., scale = TRUE, data = Penta[, -1], ncomp = 3, contr = "contr.none",
               method = "bidiagpls", validation = "oob")
smc(mod1)
plot(smc(mod1))
}