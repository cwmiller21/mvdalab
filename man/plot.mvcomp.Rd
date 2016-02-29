\name{plot.mvcomp}
\alias{plot.mvcomp}
\title{
Plot of Multivariate Mean Vector Comparison
}
\description{
Plot a comparison of mean vectors drawn from two populations.
}
\usage{
\method{plot}{mvcomp}(x, Diff2Plot = c(3, 4), segments = 51, include.zero = FALSE, ...)
}
\arguments{
\item{x}{ an plot.mvcomp object. }
\item{segments}{ number of line-segments used to draw ellipse. }
\item{Diff2Plot}{ variable differences to plot. }
\item{include.zero}{ add the zero axis to the graph output. }
\item{\dots}{ additional arguments.  Currently ignored. }
}
\details{
This function provides a plot of the T2-statistic for testing the equality of two mean vectors.  This test is appropriate for testing two populations, assuming independence.  

Assumptions:

The sample for both populations is a random sample from a multivariate population.

-Both populations are independent

-Both populations are multivariate normal

-Covariance matrices are approximately equal

If the confidence ellipse does not cover c(0, 0), we reject the NULL that the differnece between mean vectors is equal to zero (at the stated alpha level).

}
\value{
This function returns a plot of the simultaneous confidence intervals for the p-variates and its corresponding confidence ellipse at the stated confidence level.
}
\author{Nelson Lee Afanador (\email{nelson.afanador@mvdalab.com})}
\references{
Johnson, R.A., Wichern, D.W. (2002) Applied Multivariate Statistical Analysis. Prentice Hall.
}
\examples{
data(College)
dat1 <- College
#Generate a 'fake' difference of 15 units
dat2 <- College + matrix(rnorm(nrow(dat1) * ncol(dat1), mean = 15), 
        nrow = nrow(dat1), ncol = ncol(dat1))

Comparison <- MVComp(dat1, dat2, level = .95)
Comparison
plot(Comparison, Diff2Plot = c(1, 2), include.zero = FALSE)
plot(Comparison, Diff2Plot = c(1, 2), include.zero = TRUE)

plot(Comparison, Diff2Plot = c(2, 3), include.zero = FALSE)
plot(Comparison, Diff2Plot = c(2, 3), include.zero = TRUE)


data(iris)
dat1b <- iris[, -5]
#Generate a 'fake' difference of .5 units
dat2b <- dat1b + matrix(rnorm(nrow(dat1b) * ncol(dat1b), mean = .5), 
          nrow = nrow(dat1b), ncol = ncol(dat1b))

Comparison2 <- MVComp(dat1b, dat2b, level = .90)
plot(Comparison2, Diff2Plot = c(1, 2), include.zero = FALSE)
plot(Comparison2, Diff2Plot = c(1, 2), include.zero = TRUE)

plot(Comparison2, Diff2Plot = c(3, 4), include.zero = FALSE)
plot(Comparison2, Diff2Plot = c(3, 4), include.zero = TRUE)
}