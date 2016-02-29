\name{predict.mvdareg}
\alias{predict.mvdareg}
\title{Model Predictions From a plsFit Model}
\description{\code{predict} provides predictions from the results of a pls model.}
\usage{
\method{predict}{mvdareg}(object, newdata, ncomp = object$ncomp, 
            na.action = na.pass, ...) 
}
\arguments{
\item{object}{ A \code{plsFit} model. }
\item{newdata}{ An optional data frame in which to look for variables with which to predict. If omitted, the fitted values are used. }
\item{ncomp}{ the number of components to include in the model (see below). }
\item{na.action}{ function determining what should be done with missing values in newdata. The default is to predict \code{NA}. }
\item{\dots}{ additional arguments.  Currently ignored. }
}
\details{
\code{predict.mvdareg} produces predicted values, obtained by evaluating the regression function in the frame newdata (which defaults to model.frame(object).  If \code{newdata} is omitted the predictions are based on the data used for the fit.

If \code{comps} is missing (or is NULL), predictions of the number of latent variables is provided.   Otherwise, if \code{comps} is given parameters for a model with only the requested components is returned. The generic function \code{residuals} return the model residuals for all the components specified for the model. If the model was fitted with na.action = na.exclude (or after setting the default na.action to \code{na.exclude} with options), the residuals corresponding to excluded observations are returned as NA; otherwise, they are omitted.

}
\value{
\code{predict.mvdareg} produces a vector of predictions or a matrix of predictions
}
\author{Nelson Lee Afanador (\email{nelson.afanador@mvdalab.com})}
\references{
NOTE: This function is adapted from \code{mvr} in package \pkg{pls} with extensive modifications by Nelson Lee Afanador.
}
\seealso{\code{\link{coef}}, \code{\link{coefficients.boots}}, \code{\link{coefficients}}, 
\code{\link{loadings}}, \code{\link{loadings.boots}}, \code{\link{weights}}, 
\code{\link{weight.boots}}}
\examples{
mod1 <- plsFit(log.RAI ~., scale = TRUE, data = Penta[, -1], 
        ncomp = 3, contr = "contr.none", method = "bidiagpls", 
        validation = "oob")
predict.mvdareg(mod1)
residuals(mod1)
}