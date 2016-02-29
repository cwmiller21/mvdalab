\name{imputeQs}
\alias{imputeQs}
\title{Quartile Naive Imputation of Missing Values }
\description{Missing value imputed as 'Missing'.}
\usage{
imputeQs(data)
}
\arguments{
\item{data}{ a dataset with missing values }
}
\details{
A completed data frame is returned.  For continous variables with missing values, missing values are replaced with 'Missing', while the non-missing values are replaced with their corresponding quartile assignment.  For categorical variable with missing values, missing values are replaced with 'Missing'.  This procedure can greatly increases the dimensionality of the data.
}
\author{Nelson Lee Afanador (\email{nelson.afanador@mvdalab.com})}
\examples{
dat <- introNAs(iris, percent = 25)
imputeQs(dat)
}
