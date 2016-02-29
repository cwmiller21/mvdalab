\name{imputeBasic}
\alias{imputeBasic}
\title{Naive imputation of missing values.}
\description{Imputes the mean or median for continous variables; highest frequency for categorical variables.}
\usage{
imputeBasic(data, Init = "mean")
}
\arguments{
\item{data}{ a dataset with missing values }
\item{Init}{ For continous variables impute either the mean or median }
}
\details{
A completed data frame is returned. For numeric variables, \code{NAs} are replaced with column means or medians. For categorical variables, \code{NAs} are replaced with the most frequent levels. If object contains no \code{NAs}, it is returned unaltered.
}
\value{
\code{imputeBasic} returns a list containing the following components:
\item{Imputed.DataFrame}{ Final imputed data frame }
\item{Imputed.Missing.Continous}{ Imputed continous values }
\item{Imputed.Missing.Factors}{ Imputed categorical values }
}
\author{Nelson Lee Afanador (\email{nelson.afanador@mvdalab.com})}
\examples{
dat <- introNAs(iris, percent = 25)
imputeBasic(dat)
}