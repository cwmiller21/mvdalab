\name{imputeEM}
\alias{imputeEM}
\alias{print.empca}
\title{Expectation Maximization (EM) for imputation of missing values.}
\description{Missing values are iterarively updated via an EM algorithm.}
\usage{
imputeEM(data, impute.ncomps = 2, pca.ncomps = 2, CV = TRUE, Init = "mean",
         scale = TRUE, iters = 25, tol = .Machine$double.eps^0.25)
}
\arguments{
\item{data}{ a dataset with missing values. }
\item{impute.ncomps}{ integer corresponding to the minimum number of components to test. }
\item{pca.ncomps}{minimum number of components to use in the imputation.}
\item{CV}{ Use cross-validation in determining the optimal number of components to retain for the final imputation. }
\item{Init}{ For continous variables impute either the mean or median. }
\item{scale}{ Scale variables to unit variance. }
\item{iters}{ For continous variables impute either the mean or median. }
\item{tol}{ the threshold for assessing convergence. }
}
\details{
A completed data frame is returned that mirrors a \code{model.matrix}. \code{NAs} are replaced with convergence values as obtained via EM. If object contains no \code{NAs}, it is returned unaltered.  
}

\value{
\code{imputeEM} returns a list containing the following components:
\item{Imputed.DataFrames}{ A list of imputed data frames across \code{impute.comps} }
\item{Imputed.Continous}{ A list of imputed values, at each EM iteration, across \code{impute.comps} }
\item{CV.Results}{ Cross-validation results across \code{impute.comps} }
\item{ncomps}{ \code{impute.comps} }
}
\author{Nelson Lee Afanador (\email{nelson.afanador@mvdalab.com}), Thanh Tran (\email{thanh.tran@mvdalab.com})}
\references{
B. Walczak, D.L. Massart. Dealing with missing data, Part I. Chemom. Intell. Lab. Syst. 58 (2001); 15:27
}
\examples{
dat <- introNAs(iris, percent = 25)
imputeEM(dat)
}
