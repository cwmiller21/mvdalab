\name{mvdaloo}
\alias{mvdaloo}
\title{Leave-one-out routine for \code{mvdareg} objects}
\description{When \code{validation = 'loo'} this routine effects the leave-one-out cross-validation procedure for \code{mvdareg} objects.}
\usage{
mvdaloo(X, Y, ncomp, weights = NULL, method = c("bidiagpls", "pls1gm"), 
        scale = FALSE, boots = NULL, ...)
}
\arguments{
\item{X}{ a matrix of observations. \code{NAs} and Infs are not allowed. }
\item{Y}{ a vector. \code{NAs} and \code{Infs} are not allowed. }
\item{ncomp}{ the number of components to include in the model (see below). }
\item{weights}{currently not in use}
\item{method}{ PLS algorithm used }
\item{scale}{ scaling used }
\item{boots}{not applicable for \code{validation = 'loo'}}
\item{\dots}{ additional arguments.  Currently ignored. }
}
\details{This function should not be called directly, but through the generic function \code{plsFit} with the argument \code{validation = 'loo'}. 
}
\value{
Provides the following bootstrapped results as a list for \code{mvdareg} objects:
  
\item{cvR2}{ leave-one-out estimate of cvR2. }
\item{PRESS}{ leave-one-out estimate of prediction error sums of squares. }
\item{MSPRESS}{ leave-one-out estimate of mean squared error prediction sums of squares. }
\item{RMSPRESS}{ leave-one-out estimate of mean squared error prediction sums of squares. }
\item{in.bag}{ leave-one-out samples used for model building. }
}
\author{Nelson Lee Afanador (\email{nelson.afanador@mvdalab.com}), Thanh Tran (\email{thanh.tran@mvdalab.com})}
\references{

NOTE: This function is adapted from \code{mvr} in package \pkg{pls} with extensive modifications by Nelson Lee Afanador and Thanh Tran.

}
\seealso{
\code{\link{plsFit}}, \code{\link{mvdaboot}}
}
\examples{
data(Penta)
mod1 <- plsFit(log.RAI ~., scale = TRUE, data = Penta[, -1], 
        ncomp = 3, contr = "contr.none", method = "bidiagpls", 
        validation = "loo")

mod1$validation$cvR2
mod1$validation$PRESS
mod1$validation$MSPRESS
mod1$validation$RMSPRESS
mod1$validation$in.bag  
}