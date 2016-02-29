\name{mvdaboot}
\alias{mvdaboot}
\title{Bootstrapping routine for \code{mvdareg} objects}
\description{When \code{validation = 'oob'} this routine effects the bootstrap procedure for \code{mvdareg} objects.}
\usage{
mvdaboot(X, Y, ncomp, method = c("bidiagpls", "pls1gm"), 
        scale = FALSE, n_cores, boots, ...) 
}
\arguments{
\item{X}{ a matrix of observations. NAs and Infs are not allowed. }
\item{Y}{ a vector. NAs and Infs are not allowed. }
\item{ncomp}{ the number of components to include in the model (see below). }
\item{method}{ PLS algorithm used. }
\item{scale}{ scaling used. }
\item{n_cores}{ No. of cores to run for parallel processing.  Currently set to 2 (4 max). }
\item{boots}{No. of bootstrap samples when \code{validation = 'oob'}}
\item{\dots}{ additional arguments.  Currently ignored. }
}
\details{This function should not be called directly, but through the generic function plsFit with the argument \code{validation = 'oob'}. 
}
\value{
  Provides the following bootstrapped results as a list for \code{mvdareg} objects:
  
\item{coefficients}{ fitted values }
\item{weights}{ weights }
\item{loadings}{ loadings }
\item{ncomp}{ number of latent variables }
\item{bootstraps}{ No. of bootstraps }
\item{scores}{ scores }
\item{cvR2}{ bootstrap estimate of cvR2 }
\item{PRESS}{ bootstrap estimate of prediction error sums of squares }
\item{MSPRESS}{ bootstrap estimate of mean squared error prediction sums of squares }
\item{boot.means}{ bootstrap mean of bootstrapped parameters }
\item{RMSPRESS}{ bootstrap estimate of mean squared error prediction sums of squares }
\item{D2}{ bidiag2 matrix }
\item{iD2}{ Inverse of bidiag2 matrix }
\item{y.loadings}{ normalized y-loadings }
\item{y.loadings2}{ non-normalized y-loadings }
\item{MSPRESS.632}{ .632 corrected estimate of MSPRESS }
\item{oob.fitted}{ out-of-bag PLS fitted values }
\item{RMSPRESS.632}{ .632 corrected estimate of RMSPRESS }
\item{in.bag}{ bootstrap samples used for model building at each bootstrap }
}
\author{Nelson Lee Afanador (\email{nelson.afanador@mvdalab.com}), Thanh Tran (\email{thanh.tran@mvdalab.com})}
\references{
There are many references explaining the bootstrap and its implementation for confidence interval estimation. Among them are:

Davison, A.C. and Hinkley, D.V. (1997) Bootstrap Methods and Their Application. Cambridge University Press.

Efron, B. and Tibshirani, R. (1993) An Introduction to the Bootstrap. Chapman & Hall.

Hinkley, D.V. (1988) Bootstrap methods (with Discussion). Journal of the Royal Statistical Society, B, 50, 312:337, 355:370.

NOTE: This function is adapted from \code{mvr} in package \pkg{pls} with extensive modifications by Nelson Lee Afanador and Thanh Tran.

}
\seealso{
\code{\link{plsFit}}, \code{\link{mvdaloo}}
}

\examples{
data(Penta)
mod1 <- plsFit(log.RAI ~., scale = TRUE, data = Penta[, -1], 
               ncomp = 3, contr = "contr.none", method = "bidiagpls", 
               validation = "oob")

mod1$validation$coefficients
mod1$validation$weights
mod1$validation$loadings
mod1$validation$ncomp
mod1$validation$bootstraps
mod1$validation$scores
mod1$validation$cvR2
mod1$validation$PRESS
mod1$validation$MSPRESS
mod1$validation$boot.means
mod1$validation$RMSPRESS
mod1$validation$D2
mod1$validation$iD2
mod1$validation$y.loadings
mod1$validation$y.loadings2
mod1$validation$MSPRESS.632 
mod1$validation$oob.fitted
mod1$validation$RMSPRESS.632
mod1$validation$in.bag 
}