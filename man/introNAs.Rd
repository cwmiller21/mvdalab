\name{introNAs}
\alias{introNAs}
\title{Introduce NA's into a Dataframe}
\description{Function for testing missing value imputation algorithms}
\usage{
introNAs(data, percent = 25)
}
\arguments{
\item{data}{ a dataset without missing values. }
\item{percent}{ the percent data that should be randomly assigned as missing }
}
\details{
A completed data frame is returned with the desired percentage of missing data.  \code{NAs} are assigned at random.
}
\author{Nelson Lee Afanador (\email{nelson.afanador@mvdalab.com})}
\examples{
dat <- introNAs(iris)
dat
}
