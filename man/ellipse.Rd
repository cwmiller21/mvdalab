\name{ellipse}
\alias{ellipse.mvdalab}
\title{Ellipses, Data Ellipses, and Confidence Ellipses}
\description{This function draws econfidence ellipses for covariance and correlation matrices derived from from either a matrix or dataframe.}
\usage{
\method{ellipse}{mvdalab}(data, center = c(0, 0), radius = "chi", scale = TRUE, segments = 51, 
        level = c(.95, .99), ...)
}
\arguments{
\item{data}{ A dataframe }
\item{center}{ 2-element vector with coordinates of center of ellipse. }
\item{radius}{ Use of the Chi or F Distributions for setting the radius of the confidence ellipse}
\item{scale}{use correlation or covariance matrix}
\item{segments}{ number of line-segments used to draw ellipse. }
\item{level}{ draw elliptical contours at these (normal) probability or confidence levels. }
\item{\dots}{ additional arguments.  Currently ignored. }
}
\details{
\code{ellipse} uses the singular value decomposition in order to generate the desired confidence regions. The default confidence ellipse is based on the chisquare statistic.
}
\author{Nelson Lee Afanador (\email{nelson.afanador@mvdalab.com})}
\references{

Fox, J. (2008) Applied Regression Analysis and Generalized Linear Models, Second Edition. Sage.

Fox, J. and Weisberg, S. (2011) An R Companion to Applied Regression, Second Edition, Sage.

}
\examples{
data(iris)
ellipse.mvdalab(iris)
}