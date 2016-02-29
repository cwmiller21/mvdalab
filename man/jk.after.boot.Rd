\name{jk.after.boot}
\alias{jk.after.boot}
\title{Jackknife After Bootstrap}
\description{This function calculates the jackknife influence values from a bootstrap output \code{mvdareg} object and plots the corresponding jackknife-after-bootstrap plot.}
\usage{
jk.after.boot(object, ncomp = object$ncomp, 
              type = c("coefficients", "loadings", "weights"), 
              parm = NULL)
}
\arguments{
\item{object}{ an mvdareg object.  A fitted model.  }
\item{ncomp}{ the component number to include in the jackknife-after-bootstrap plot assessment. }
\item{type}{ input parameter vector.  }
\item{parm}{ predictor variable for which to perform the assessment.  if \code{NULL} one will be chosen at random.}
}
\details{
The centred jackknife quantiles for each observation are estimated from those bootstrap samples in which a particular observation did not appear. These are then plotted against the influence values. 

The resulting plots are useful diagnostic tools for looking at the way individual observations affect the bootstrap output.

The plot will consist of a number of horizontal dotted lines which correspond to the quantiles of the centred bootstrap distribution. For each data point the quantiles of the bootstrap distribution calculated by omitting that point are plotted against the jackknife values. The observation number is printed below the plots. To make it easier to see the effect of omitting points on quantiles, the plotted quantiles are joined by line segments. These plots provide a useful diagnostic tool in establishing the effect of individual observations on the bootstrap distribution. See the references below for some guidelines on the interpretation of the plots.
}
\value{
There is no returned value but a graph is generated on the current graphics display.
}
\author{Nelson Lee Afanador (\email{nelson.afanador@mvdalab.com})}
\references{
Davison, A.C. and Hinkley, D.V. (1997) Bootstrap Methods and Their Application. Cambridge University Press.

Efron, B. (1992) Jackknife-after-bootstrap standard errors and influence functions (with Discussion). Journal of the Royal Statistical Society, B, 54, 83:127.
}
\examples{
mod1 <- plsFit(log.RAI ~., scale = TRUE, data = Penta[, -1], 
        ncomp = 3, contr = "contr.none", method = "bidiagpls", 
        validation = "oob")
jk.after.boot(mod1, type = "loadings")
jk.after.boot(mod1, type = "weights")
jk.after.boot(mod1, type = "coefficients")
}