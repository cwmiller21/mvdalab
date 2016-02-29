\name{SeqimputeEM}
\alias{SeqimputeEM}
\alias{print.seqem}
\title{Sequential Expectation Maximization (EM) for imputation of missing values.}
\description{Missing values are sequentially updated via an EM algorithm.}
\usage{
SeqimputeEM(data, max.ncomps = 5, max.ssq = 0.99, Init = "mean", 
            adjmean = FALSE, max.iters = 200, 
            tol = .Machine$double.eps^0.25)
}
\arguments{
\item{data}{ a dataset with missing values. }
\item{max.ncomps}{ integer corresponding to the maximum number of components to test }
\item{max.ssq}{ maximal SSQ for final number of components. This will be improved by automation.  }
\item{Init}{ For continous variables impute either the mean or median. }
\item{adjmean}{ Adjust (recalculate) mean after each iteration. }
\item{max.iters}{ maximum number of iterations for the algorithm. }
\item{tol}{ the threshold for assessing convergence. }
}
\details{
A completed data frame is returned that mirrors the model matrix. \code{NAs} are replaced with convergence values as obtained via Seqential EM algorithm. If object contains no \code{NAs}, it is returned unaltered.  
}

\value{
\item{Imputed.DataFrames}{ A list of imputed data frames across \code{impute.comps} }
\item{ncomps}{ number of components to test }
}
\author{Thanh Tran (\email{thanh.tran@mvdalab.com}), Nelson Lee Afanador (\email{nelson.afanador@mvdalab.com})}
\references{
NOTE: Publication Pending
}
\examples{
dat <- introNAs(iris, percent = 25)
SeqimputeEM(dat)
}