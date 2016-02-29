\name{smc.acfTest}
\alias{smc.acfTest}
\title{Test of the Residual Significant Multivariate Correlation Matrix for the presence of Autocorrelation}
\description{This function peforms a 1st order test of the Residual Significant Multivariate Correlation Matrix in order to help determine if the \code{smc} should be performed correcting for 1st order autocorrelation.}
\usage{
smc.acfTest(object, ncomp = object$ncomp)
}
\arguments{
\item{object}{ an object of class \code{mvdareg}, i.e. \code{plsFit}. }
\item{ncomp}{ the number of components to include in the acf assessment }
}
\details{
This function computes a test for 1st order auto correlation in the \code{smc} residual matrix.
}
\value{
The output of \code{smc.acfTest} is a list detailing the following:

\item{variable}{ variable for whom the test is being performed }
\item{ACF}{ value of the 1st lag of the ACF }
\item{Significant}{ Assessment of the statistical significance of the 1st order lag }

}
\author{Nelson Lee Afanador (\email{nelson.afanador@mvdalab.com})}
\references{
Thanh N. Tran, Nelson Lee Afanador, Lutgarde M.C. Buydens, Lionel Blanchet, Interpretation of variable importance in Partial Least Squares with Significance Multivariate Correlation (sMC). Chemom. Intell. Lab. Syst. 2014; 138: 153:160.

Nelson Lee Afanador, Thanh N. Tran, Lionel Blanchet, Lutgarde M.C. Buydens, Variable importance in PLS in the presence of autocorrelated data - Case studies in manufacturing processes. Chemom. Intell. Lab. Syst. 2014; 139: 139:145.
}
\examples{
data(Penta)
mod1 <- plsFit(log.RAI ~., scale = TRUE, data = Penta[, -1], ncomp = 3, contr = "contr.none",
               method = "bidiagpls", validation = "oob")
smc.acfTest(mod1, ncomp = 3)
}