\name{MVComp}
\alias{MVComp}
\alias{print.mvcomp}
\title{
Traditional Multivariate Mean Vector Comparison
}
\description{
Performs a traditional multivariate comparison of mean vectors drawn from two populations.
}
\usage{
MVComp(data1, data2, level = .95)
}
\arguments{
\item{data1}{ a multivariable dataset to compare to. }
\item{data2}{ a multivariable dataset to compare. }
\item{level}{ draw elliptical contours at these (normal) probability or confidence levels. }
}
\details{
This function provides a T2-statistic for testing the equality of two mean vectors.  This test is appropriate for testing two populations, assuming independence.  

Assumptions:

The sample for both populations is a random sample from a multivariate population.

-Both populations are independent

-Both populations are multivariate normal

-Covariance matrices are approximately equal

}
\value{
This function returns the simultaneous confidence intervals for the p-variates and its corresponding confidence ellipse at the stated confidence level.
}
\references{
Johnson, R.A., Wichern, D.W. (2002) Applied Multivariate Statistical Analysis. Prentice Hall.
}
\author{Nelson Lee Afanador (\email{nelson.afanador@mvdalab.com})}
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