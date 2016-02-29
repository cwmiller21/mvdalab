\name{vip}
\alias{vip}
\alias{vip.boots}
\alias{print.vip}
\title{Variable Importance in the Projection}
\description{This function calculated the variable importance in the projection (VIP) metric for an \code{mvdareg} object}
\usage{
vip(object, ncomp = object$ncomp, conf = .95)
}
\arguments{
\item{object}{ an mvdareg or mvdapaca object.  A fitted model. }
\item{ncomp}{ the number of components to include in the model (see below). }
\item{conf}{ for a bootstrapped model, the confidence level to use. }
}
\details{
\code{vip} is used to extract a summary of the variable importance in the projection of a PLS model.

If \code{comps} is missing (or is \code{NULL}), summaries for all regression estimates are returned.  Otherwise, if \code{comps} are given parameters for a model with only the requested component comps is returned. 

For \code{mvdareg} objects only, boostrap summaries provided are for actual VIPs, bootstrap percentiles, bootstrap mean, skewness, and bias.  
}
\value{
A \code{vip} object contains component results for the following:
\item{ncomp}{ the number of components to include in the model. }
\item{variable}{ variable names. }
\item{actual}{ Actual loading estimate using all the data. }
\item{percentiles}{ confidence intervals. }
\item{boot.mean}{ mean of the bootstrap. }
\item{skewness}{ skewness of the bootstrap distribution. }
\item{bias}{ estimate of bias w.r.t. the loading estimate. }
\item{Bootstrap Error}{ estimate of bootstrap standard error. }
\item{t value}{ approximate 't-value' based on the \code{Bootstrap Error}. }
\item{bias corrected t value}{ approximate 'bias corrected t-value' based on the \code{Bootstrap Error}. }
}
\author{Nelson Lee Afanador (\email{nelson.afanador@mvdalab.com})}
\references{
Il-Gyo Chong, Chi-Hyuck Jun, Performance of some variable selection methods when multicollinearity. Chemom. Intell. Lab. Syst. 2004; 78: 103:112.

Nelson Lee Afanador, Thanh N. Tran, Lutgarde M.C. Buydens, An assessment of the jackknife and bootstrap procedures on uncertainty estimation in the variable importance in the projection metric. Chemom. Intell. Lab. Syst. 2014; 137: 162:172.
}
\seealso{\code{\link{smc}}, \code{\link{sr}}}
\examples{
data(Penta)
mod1 <- plsFit(log.RAI ~., scale = TRUE, data = Penta[, -1], ncomp = 3, contr = "contr.none",
               method = "bidiagpls", validation = "oob")
vip(mod1, conf = .95)
vip(mod1, ncomp = 2, conf = .95)
vip(mod1, ncomp = 3, conf = .95)


mod2 <- plsFit(log.RAI ~., scale = TRUE, data = Penta[, -1], ncomp = 3, contr = "contr.none",
               method = "bidiagpls", validation = "loo")
vip(mod2, ncomp = 1, conf = .95)
vip(mod2, ncomp = 2, conf = .95)
vip(mod2, ncomp = 3, conf = .95)
}