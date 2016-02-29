\name{my.dummy.df}
\alias{my.dummy.df}
\title{Create a Design Matrix with the Desired Constrasts}
\description{This function generates a dummy variable data frame in support various functions.}
\usage{my.dummy.df(data, contr = "contr.none")}
\arguments{
\item{data}{ a data frame }
\item{contr}{ an optional list. See the contrasts.arg of model.matrix.default. }
}
\details{
my.dummy.df takes a \code{data.frame} with categorical variables, and returns a \code{data.frame} in which all the categorical variables columns are expanded as dummy variables. 

The argument \code{contr} is passed to the default \code{contr.none}; \code{contr.helmert}, \code{contr.poly}, \code{contr.sum}, \code{contr.treatment} are also supported.
}
\value{
For datasets with categorical variables it produces the specified design matrix.
}
\author{Nelson Lee Afanador (\email{nelson.afanador@mvdalab.com})}
\examples{
data(iris)
my.dummy.df(iris)
}