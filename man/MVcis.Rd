\name{MVcis}
\alias{MVcis}
\title{
Calculate Hotelling's T2 Confidence Intervals
}
\description{
Calculate joint confidence intervals (Hotelling's T2 Intervals).
}
\usage{
MVcis(data, segments = 51, level = .95, Vars2Plot = c(1, 2), include.zero = F)
}
\arguments{
\item{data}{ a multivariable dataset to compare to means}
\item{segments}{ number of line-segments used to draw ellipse. }
\item{level}{ draw elliptical contours at these (normal) probability or confidence levels. }
\item{Vars2Plot}{ variables to plot }
\item{include.zero}{ add the zero axis to the graph output }
}
\details{
This function calculates the Hotelling's T2 Intervals for a mean vector.

Assumption:

Population is a random sample from a multivariate population.

If the confidence ellipse does not cover c(0, 0), we reject the NULL that the joint confidence region is equal to zero (at the stated alpha level).

}
\value{
This function returns the Hotelling's T2 confidence intervals for the p-variates and its corresponding confidence ellipse at the stated confidence level.
}
\references{
Johnson, R.A., Wichern, D.W. (2002) Applied Multivariate Statistical Analysis. Prentice Hall.
}
\author{Nelson Lee Afanador (\email{nelson.afanador@mvdalab.com})}
\seealso{
\code{\link{MVComp}}
}
\examples{
data(College)
MVcis(College, Vars2Plot = c(1, 2), include.zero = TRUE)
}