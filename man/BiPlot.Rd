\name{BiPlot}
\alias{BiPlot}
\title{Generates a biplot from the output of an 'mvdareg' and 'mvdapca' object}
\description{Generates a 2D Graph of both the scores and loadings for both \code{"mvdareg"} and \code{"mvdapca"} objects.}
\usage{
BiPlot(object, diag.adj = c(0, 0), axis.scaling = 2, cov.scale = FALSE)
}
\arguments{
\item{object}{ an object of class \code{"mvdareg"} or \code{"mvdapca"}. }
\item{diag.adj}{ adjustment to singular values.  see details. }
\item{axis.scaling}{ a graphing parameter for extenting the axis. }
\item{cov.scale}{implement covariance scaling}
}
\details{
\code{"BiPlot"} is used to extract a 2D graphical summary of the scores and loadings of PLS and PCA models.  

The singular values are scaled so that the approximation becomes X = GH':

X = ULV' = (UL^alpha1)(L^alpha2V') = GH', and where alpha2 is = to (1 = alpha)

The rows of the G matrix are plotted as points, corresponding to observations. The rows of the H matrix are plotted as vectors, corresponding to variables.  The choice of alpha determines the following:

c(0, 0):  variables are scaled to unit length and treats observations and variables symmetrically. 

c(0, 1):  This biplot attempts to preserve relationships between variables wherein the distance betweein any two rows of G is proportional to the Mahalanobis distance between the same observations in the data set.
 
c(1, 0):  This biplot attempts to preserve the distance between observations where in the positions of the points in the biplot are identical to the score plot of first two principal components, but the distance between any two rows of G is equal to the Euclidean distance between the corresponding observations in the data set.

\code{cov.scale = FALSE} sets diag.adj to c(0, 0) and multiples G by sqrt(n - 1) and divides H by sqrt(n - 1).  In this biplot the rows of H approximate the variance of the corresponding variable, and the distance between any two points of G approximates the Mahalanobis distance between any two rows.

Additional scalings may be implemented.

}
\author{Nelson Lee Afanador (\email{nelson.afanador@mvdalab.com})}

\references{

SAS Stat Studio 3.11 (2009), User's Guide.

Additional information pertaining to biplots can be obtained from the following:

Friendly, M. (1991), SAS System for Statistical Graphics , SAS Series in Statistical Applications, Cary, NC: SAS Institute

Gabriel, K. R. (1971), "The Biplot Graphical Display of Matrices with Applications to Principal Component Analysis," Biometrika , 58(3), 453--467. 

Golub, G. H. and Van Loan, C. F. (1989), Matrix Computations , Second Edition, Baltimore: Johns Hopkins University Press. 

Gower, J. C. and Hand, D. J. (1996), Biplots , London: Chapman & Hall. 

Jackson, J. E. (1991), A User's Guide to Principal Components , New York: John Wiley & Sons.

}

\examples{
data(Penta)
mod1 <- plsFit(log.RAI ~., scale = TRUE, data = Penta[, -1], 
               ncomp = 3, contr = "contr.none",
               method = "bidiagpls", validation = "oob")
BiPlot(mod1, diag.adj = c(0, 0), axis.scaling = 2, cov.scale = FALSE)
BiPlot(mod1, diag.adj = c(1, 0), axis.scaling = 2, cov.scale = FALSE)
BiPlot(mod1, diag.adj = c(0, 1), axis.scaling = 2, cov.scale = FALSE)
BiPlot(mod1, axis.scaling = 2, cov.scale = TRUE)

data(Penta)
mod2 <- pcaFit(Penta[, -1])
BiPlot(mod2, diag.adj = c(0, 0), axis.scaling = 2.25, cov.scale = FALSE)
BiPlot(mod2, diag.adj = c(1, 0), axis.scaling = 2.25, cov.scale = FALSE)
BiPlot(mod2, diag.adj = c(0, 1), axis.scaling = 2.25, cov.scale = FALSE)
BiPlot(mod2, axis.scaling = 2.25, cov.scale = TRUE)
}