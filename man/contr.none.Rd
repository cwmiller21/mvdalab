\name{contr.none}
\alias{contr.none}
\title{Cell Means Contrast Matrix }
\description{This function generates a cell means contrast matrix to support various functions.}
\usage{contr.none(n, contrasts)}
\arguments{
\item{n}{ A vector of levels for a factor, or the number of levels. }
\item{contrasts}{ A logical indicating whether contrasts should be computed. This argument is ignored in contr.none. }
}
\details{
This function,as authored by Jelle Goeman, has been imported from the \pkg{penalized} package and generates a cell means contrast matrix in support of various functions. 
}
\value{
For datasets with categorical variables it produces the needed design matrix.
}
\author{Jelle Goeman}
\references{
Original:  This very useful function was obtained from the \pkg{penalized} package and has been imported to prevent additional loading time.  Full credit and thanks are given to the original author, Jelle Goeman.
}
\examples{
# Three levels
levels <- LETTERS[1:3]
contr.none(levels)

# Two levels
levels <- LETTERS[1:2]
contr.none(levels)
}