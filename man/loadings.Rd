\name{loadings}
\alias{loadings.mvdareg}
\title{Summary Information Pertaining to the Bootstrapped Loadings}
\description{Functions to extract loadings bootstrap information from mvdalab objects.}
\usage{
\method{loadings}{mvdareg}(object, ncomp = object$ncomp, conf = .95, ...)
}
\arguments{
\item{object}{ an mvdareg or mvdapaca object.  A fitted model. }
\item{ncomp}{ the number of components to include in the model (see below). }
\item{conf}{ for a bootstrapped model, the confidence level to use. }
\item{\dots}{ additional arguments.  Currently ignored. }
}
\details{
\code{loadings} is used to extract a summary of the loadings of a PLS or PCA model.
If \code{ncomps} is missing (or is NULL), summaries for all loadings estimates are returned.  Otherwise, if comps is given parameters for a model with only the requested component comps is returned. 

Boostrap summaries are provided for \code{mvdareg} objects where \code{validation = "oob"}.  These summaries can also be extracted using \code{loadings.boots}
}
\value{
A loadings object contains a data frame with columns:
\item{variable}{ variable names }
\item{Actual}{ Actual loading estimate using all the data }
\item{BCa percentiles}{ confidence intervals }
\item{boot.mean}{ mean of the bootstrap }
\item{skewness}{ skewness of the bootstrap distribution }
\item{bias}{ estimate of bias w.r.t. the loading estimate }
\item{Bootstrap Error}{ estimate of bootstrap standard error }
\item{t value}{ approximate 't-value' based on the \code{Bootstrap Error} }
\item{bias t value}{ approximate 'bias t-value' based on the \code{Bootstrap Error} }
}
\author{Nelson Lee Afanador (\email{nelson.afanador@mvdalab.com})}
\seealso{\code{\link{loadingsplot}}, \code{\link{loadings.boots}}, \code{\link{loadingsplot2D}}}
\references{
There are many references explaining the bootstrap. Among them are:

Davison, A.C. and Hinkley, D.V. (1997) Bootstrap Methods and Their Application. Cambridge University Press.

Efron, B. (1992) Jackknife-after-bootstrap standard errors and influence functions (with Discussion). Journal of the Royal Statistical Society, B, 54, 83:127.
}
\examples{
data(Penta)
mod1 <- plsFit(log.RAI ~., scale = TRUE, data = Penta[, -1], 
               ncomp = 3, contr = "contr.none", method = "bidiagpls", 
               validation = "oob")
loadings(mod1, ncomp = 3, conf = .95)

data(iris)
pc1 <- pcaFit(iris)
loadings(pc1)
}