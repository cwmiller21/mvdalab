\name{imputeRough}
\alias{imputeRough}
\alias{print.roughImputation}
\title{Naive Imputation of Missing Values for Dummy Variable Model Matrix}
\description{After generating a cell means model matrix, impute expected values (mean or median for continous; hightest frequency for categorical).}
\usage{imputeRough(data, Init = "mean")}
\arguments{
\item{data}{ a dataset with missing values }
\item{Init}{ For continous variables impute either the mean or median }
}
\details{
A completed data frame is returned that mirrors a \code{model.matrix}. \code{NAs} are replaced with column means or medians. If object contains no \code{NAs}, it is returned unaltered.  This is the starting point for imputeEM.
}
\value{
\code{imputeRough} returns a list containing the following components:

\item{Initials}{ Imputed values }
\item{Pre.Imputed}{ Pre-imputed data frame }
\item{Imputed.Dataframe}{ Imputed data frame }
}
\author{Nelson Lee Afanador (\email{nelson.afanador@mvdalab.com})}
\examples{
dat <- introNAs(iris, percent = 25)
imputeRough(dat)
}
