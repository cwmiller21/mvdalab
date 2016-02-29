\name{coefficients.mvdareg}
\alias{coefficients.mvdareg}
\title{Extract Summary Information Pertaining to the Coefficients resulting from a PLS model}
\description{Functions to extract regression coefficient bootstrap information from mvdalab objects.}
\usage{\method{coefficients}{mvdareg}(object, ncomp = object$ncomp, conf = .95, ...)}
\arguments{
\item{object}{ an mvdareg object.  A fitted model. }
\item{ncomp}{ the number of components to include in the model (see below). }
\item{conf}{ for a bootstrapped model, the confidence level to use. }
\item{\dots}{ additional arguments.  Currently ignored. }
}
\details{

\code{coefficients} is used to extract a bootstrap summary of the regression of a PLS model.

If \code{comps} is missing (or is NULL), summaries for all regression estimates are returned.  Otherwise, if comps is given parameters for a model with only the requested component comps is returned. 

Boostrap summaries provided are for actual regression coefficients, bootstrap percentiles, bootstrap mean, skewness, and bias.  These summaries can also be extracted using \code{coefficients.boots}

}

\value{
A coefficients object contains a data frame with columns:
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
\seealso{\code{\link{coef}}, \code{\link{coefficients.boots}}, \code{\link{coefficients}}}
\examples{
data(Penta)
mod1 <- plsFit(log.RAI ~., scale = TRUE, data = Penta[, -1], 
        ncomp = 3, contr = "contr.none", method = "bidiagpls", 
        validation = "oob")
coefficients(mod1)
}