\name{delete.intercept}
\alias{delete.intercept}
\title{Delete Intercept from Model Matrix}
\description{A utility function from the \pkg{pls} package to delete any intercept column from a model matrix, and adjust the "assign" attribute correspondingly.}
\usage{delete.intercept(mm)}
\arguments{\item{mm}{ Model Matrix }}
\value{
A model matrix without intercept column.
}
\author{Bjorn-Helge Mevik and Ron Wehrens}
\references{
Original: This very useful function was obtained from the \pkg{pls} package and has been imported to prevent additional loading time. Full credit and thanks are given to the original author, Bjorn-Helge Mevik and Ron Wehrens.
}