\name{plsFit}
\alias{plsFit}
\alias{mvdareg}
\alias{summary.mvdareg}
\alias{summary.mvdareg.default}
\title{Partial Least Squares Regression}
\description{Functions to perform partial least squares regression with a formula interface. Bootstraping can be used. Prediction, residuals, model extraction, plot, print and summary methods are also implemented.}
\usage{
plsFit(formula, ncomp, data, subset, na.action, contr = "contr.none",
        method = c("bidiagpls", "pls1gm"), scale = TRUE, n_cores = 2, 
        validation = c("none", "oob", "loo"), boots = 1000, model = TRUE, 
        x = FALSE, y = FALSE, ...)
  
\method{summary}{mvdareg}(object, ncomp = object$ncomp, digits = 3, ...)
}
\arguments{
\item{formula}{ a model formula (see below). }
\item{ncomp}{ the number of components to include in the model (see below). }
\item{data}{ an optional data frame containing the variables in the model. }
\item{subset}{ an optional vector specifying a subset of observations to be used in the fitting process. }
\item{na.action}{ a function which indicates what should happen when the data contain \code{NAs}. The default is set by the \code{na.action} setting of options, and is \code{na.fail} if that is unset. The default is \code{na.omit}. Another possible value is \code{NULL}, no action. Value \code{na.exclude} can be useful. }
\item{contr}{ an optional list. See the \code{contrasts.arg} of \code{model.matrix.default}. }
\item{method}{ the multivariate regression algorithm to be used. }
\item{scale}{ an optional data frame containing the variables in the model. }
\item{n_cores}{ Number of cores to run for parallel processing.  Currently set to 2 with the max being 4. }
\item{validation}{ character. What kind of (internal) validation to use. See below. }
\item{boots}{Number of bootstrap samples when \code{validation = 'oob'}}
\item{model}{ an optional data frame containing the variables in the model. }
\item{x}{ a logical. If TRUE, the model matrix is returned. }
\item{y}{ a logical. If TRUE, the response is returned. }
\item{object}{ an object of class \code{"mvdareg"}, i.e., a fitted model. }
\item{digits}{the number of decimal place to output with \code{summary.mvdareg}}
\item{\dots}{ additional arguments, passed to the underlying fit functions, and \code{mvdareg}. Currently not in use.}
}
\details{The function fits a partial least squares (PLS) model with 1, ..., \code{ncomp} number of latent variables. Multi-response models are not supported.

The type of model to fit is specified with the method argument. Two PLS algorithms are available: the bigiag2 algorithm ("bigiagpls") and the Gram-Schmidt classical orthogonal scores algorithm ("pls1gm"). 

The formula argument should be a symbolic formula of the form response ~ terms, where response is the name of the response vector and terms is the name of one or more predictor matrices, usually separated by +, e.g., y ~ X + Z. See \code{\link{lm}} for a detailed description. The named variables should exist in the supplied data data frame or in the global environment. The chapter Statistical models in R of the manual An Introduction to R distributed with R is a good reference on formulas in R.

The number of components to fit is specified with the argument \code{ncomp}. It this is not supplied, the maximal number of components is used.

If \code{validation = "oob"}, bootstrap cross-validation is performed. Bootstrap confidence intervals are provided for \code{\link{coefficients}}, \code{weights}, \code{loadings}, and \code{y.loadings}.  The number of bootstrap samples is specified with the argument \code{boots}. See \code{mvdaboot} for details. If \code{validation = "loo"}, leave-one-out cross-validation is performed. If \code{validation = "none"}, no cross-validation is performed.

The argument \code{contr} is passed to the default \code{contr.none}; \code{contr.helmert}, \code{contr.poly}, \code{contr.sum}, \code{contr.treatment} are also supported.
}
\value{
An object of class \code{mvdareg} is returned. The object contains all components returned by the underlying fit function. In addition, it contains the following:
\item{loadings}{ X loadings }
\item{weights}{ weights }
\item{D2.values}{ bidiag2 matrix }
\item{iD2}{ inverse of bidiag2 matrix }
\item{Ymean}{ mean of reponse variable }
\item{Xmeans}{ mean of predictor variables }
\item{coefficients}{ PLS regression coefficients }
\item{y.loadings}{ y-loadings }
\item{scores}{ X scores }
\item{R}{ orthogonal weights }
\item{Y.values}{ scaled response values }
\item{Yactual}{ actual response values }
\item{fitted}{ fitted values }
\item{residuals}{ residuals }
\item{Xdata}{ X matrix }
\item{iPreds}{ predicted values }
\item{y.loadings2}{ scaled y-loadings }
\item{ncomp}{ number of latent variables }
\item{method}{ PLS algorithm used }
\item{scale}{ scaling used }
\item{validation}{ validation method }
\item{call}{ model call }
\item{terms}{ model terms }
\item{model}{ fitted model }
}
\author{Nelson Lee Afanador (\email{nelson.afanador@mvdalab.com}), Thanh Tran (\email{thanh.tran@mvdalab.com})}
\references{
NOTE: This function is adapted from \code{mvr} in package \pkg{pls} with extensive modifications by Nelson Lee Afanador and Thanh Tran.
}
\seealso{\code{\link{bidiagpls.fit}}, \code{\link{mvdaboot}}, \code{\link{boot.plots}},
\code{\link{R2s}}, \code{\link{PE}}, \code{\link{ap.plot}}, 
\code{\link{T2}}, \code{\link{Xresids}}, \code{\link{smc}}, 
\code{\link{scoresplot}}, \code{\link{ScoreContrib}}, \code{\link{sr}}, 
\code{\link{loadingsplot}}, \code{\link{weightsplot}}, \code{\link{coefsplot}}, 
\code{\link{loadingsplot2D}}, \code{\link{weightsplot2D}}, \code{\link{vip}}, 
\code{\link{bca.cis}}, \code{\link{coefficients.boots}}, \code{\link{loadings.boots}}, 
\code{\link{weight.boots}}, \code{\link{coefficients}}, \code{\link{loadings}}, 
\code{\link{weights}}, \code{\link{BiPlot}}, \code{\link{jk.after.boot}}}
\examples{
###  PLS MODEL FIT WITH validation = 'oob', i.e. bootstrapping ###
data(Penta)
mod1 <- plsFit(log.RAI ~., scale = TRUE, data = Penta[, -1], 
               ncomp = 3, contr = "contr.none", method = "bidiagpls", 
               validation = "oob")

summary(mod1) #Model summary

R2s(mod1) #R2's

plot(R2s(mod1)) #R2's plot

PE(mod1) #X-explained variance

ap.plot(mod1, ncomp = 1) #actual vs. predicted plot for 1 LV
ap.plot(mod1, ncomp = 2) #actual vs. predicted plot for 2 LV
ap.plot(mod1, ncomp = 3) #actual vs. predicted plot for 3 LV

predict(mod1, ncomp = 1:3)
residuals(mod1)

loadings.boots(mod1)

boot.plots(mod1, type = "coefs", parm = NULL)
boot.plots(mod1, type = "weights", parm = NULL)
boot.plots(mod1, type = "loadings", parm = NULL)

bca.cis(mod1, conf = .95, type = "coefficients")
bca.cis(mod1, conf = .95, type = "loadings")
bca.cis(mod1, conf = .95, type = "weights")

loadingsplot(mod1, ncomp = 1, conf = 0.95) #loadings plot
weightsplot(mod1, ncomp = 2, conf = 0.95) #weights plot
coefsplot(mod1, ncomp = 3, conf = 0.95) #coef plot

coefficients(mod1, ncomp = 1, conf = .95)
loadings(mod1, ncomp = 1:2, conf = .95)
weights(mod1, ncomp = 3, conf = .95)
y.loadings(mod1, conf = .95)

jk.after.boot(mod1, type = "loadings", parm = NULL)
jk.after.boot(mod1, type = "weights", parm = NULL)
jk.after.boot(mod1, type = "coefficients", parm = NULL)

T2(mod1, ncomp = 2) #T2 plot

Xresids(mod1, ncomp = 2) #X-residuals plot
XresidualContrib(mod1, obs1 = 1)

(SMC <- smc(mod1, ncomp = 2, corrected = FALSE)) #smc variable importance
plot(SMC) #smc variable importance plot

(VIP <- vip(mod1, ncomp = 3)) #VIP variable importance
plot(VIP, ncomp = 1:3) #VIP variable importance plot

(SR <- sr(mod1, ncomp = 2)) #Selectivity ratio variable importance
plot(SR) #Plot Selectivity Ratio variable importance

scoresplot(mod1) #scoresplot variable importance

(SC <- ScoreContrib(mod1, obs1 = 1:9, obs2 = 10:11)) #score contribution
plot(SC)  #score contribution plot

loadingsplot2D(mod1, comps = c(1, 2)) #2-D loadings plot
loadingsplot2D(mod1, comps = c(2, 3)) #2-D loadings plot

weightsplot2D(mod1, comps = c(1, 2)) #2-D weights plot
weightsplot2D(mod1, comps = c(2, 3)) #2-D weights plot

BiPlot(mod1, diag.adj = c(0, 0), axis.scaling = 2, cov.scale = FALSE)
BiPlot(mod1, diag.adj = c(1, 0), axis.scaling = 2, cov.scale = FALSE)
BiPlot(mod1, diag.adj = c(0, 1), axis.scaling = 2, cov.scale = FALSE)
BiPlot(mod1, axis.scaling = 2, cov.scale = TRUE)

###  PLS MODEL FIT WITH validation = 'loo', i.e. leave-one-out CV ###

mod2 <- plsFit(log.RAI ~., scale = TRUE, data = Penta[, -1], 
               ncomp = 3, contr = "contr.none", method = "bidiagpls", 
               validation = "loo")
               
summary(mod2) #Model summary

R2s(mod2) #R2's

plot(R2s(mod2)) #R2's plot

PE(mod2) #X-explained variance

loadingsplot(mod2, ncomp = 1, conf = 0.95) #loadings plot
weightsplot(mod2, ncomp = 2, conf = 0.95) #weights plot
coefsplot(mod2, ncomp = 3, conf = 0.95) #coef plot

coefficients(mod2, ncomp = 1, conf = .95)
loadings(mod2, ncomp = 1:2, conf = .95)
weights(mod2, ncomp = 3, conf = .95)
y.loadings(mod2, conf = .95)

ap.plot(mod2, ncomp = 1) #actual vs. predicted plot for 1 LV
ap.plot(mod2, ncomp = 2) #actual vs. predicted plot for 2 LV
ap.plot(mod2, ncomp = 3) #actual vs. predicted plot for 3 LV

predict(mod2, ncomp = 1:3)
residuals(mod2)

T2(mod2, ncomp = 2) #T2 plot

Xresids(mod2, ncomp = 2) #X-residuals plot
XresidualContrib(mod2, obs1 = 1)

(SMC <- smc(mod2, ncomp = 2, corrected = FALSE)) #smc variable importance
plot(SMC) #smc variable importance plot

(VIP <- vip(mod2, ncomp = 3)) #VIP variable importance
plot(VIP, ncomp = 1:3) #VIP variable importance plot

(SR <- sr(mod2, ncomp = 2)) #Selectivity ratio variable importance
plot(SR) #Plot Selectivity Ratio variable importance

scoresplot(mod2) #scoresplot variable importance

(SC <- ScoreContrib(mod2, obs1 = 1:9, obs2 = 10:11))  #score contribution
plot(SC)  #score contribution plot

loadingsplot2D(mod2, comps = c(1, 2)) #2-D loadings plot
loadingsplot2D(mod2, comps = c(2, 3)) #2-D loadings plot

weightsplot2D(mod2, comps = c(1, 2)) #2-D weights plot
weightsplot2D(mod2, comps = c(2, 3)) #2-D weights plot

BiPlot(mod2, diag.adj = c(0, 0), axis.scaling = 2, cov.scale = FALSE)
BiPlot(mod2, diag.adj = c(1, 0), axis.scaling = 2, cov.scale = FALSE)
BiPlot(mod2, diag.adj = c(0, 1), axis.scaling = 2, cov.scale = FALSE)
BiPlot(mod2, axis.scaling = 2, cov.scale = TRUE)
}