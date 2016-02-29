\name{ScoreContrib}
\alias{ScoreContrib}
\title{Generates a score contribution plot}
\description{Generates a the Score Contribution Graph both \code{mvdareg} and \code{mvdapca} objects.}
\usage{ScoreContrib(object, obs1 = 1, obs2 = NULL)}
\arguments{
\item{object}{ an object of class \code{mvdareg} or \code{mvdapca}. }
\item{obs1}{ the first observaion(s) in the score(s) comparison. }
\item{obs2}{ the second observaion(s) in the score(s) comparison. }

}
\details{
\code{ScoreContrib} is used to generates the score contributions for both PLS and PCA models.  Up to two groups of score(s) can be selected.  If only one group is selected, the contribution is measured to the model average.  For PLS models the PCA loadings are replaced with the PLS weights.
}
\value{
The output of \code{ScoreContrib} is a matrix of score contributions for the specified observation(s).
}
\author{Nelson Lee Afanador (\email{nelson.afanador@mvdalab.com})}
\references{
MacGregor, Process Monitoring and Diagnosis by Multiblock PLS Methods, May 1994 Vol. 40, No. 5 AIChE Journal.
}
\examples{
#PLS Model
mod1 <- plsFit(log.RAI ~., scale = TRUE, data = Penta[, -1], ncomp = 3, contr = "contr.none",
               method = "bidiagpls", validation = "oob")
Score.Contributions1 <- ScoreContrib(mod1, obs1 = 1, obs2 = 3)
plot(Score.Contributions1)
Score.Contributions2 <- ScoreContrib(mod1, obs1 = c(1, 3), obs2 = c(5:10))
plot(Score.Contributions2)
Score.Contributions3 <- ScoreContrib(mod1, obs1 = 1:10, obs2 = 11:15)
plot(Score.Contributions3)
Score.Contributions4 <- ScoreContrib(mod1, obs1 = 1, obs2 = 3)
plot(Score.Contributions4)


#PCA Model
pc1 <- pcaFit(Penta[, -1])
Score.Contributions1 <- ScoreContrib(pc1, obs1 = 1, obs2 = 3)
plot(Score.Contributions1)
Score.Contributions2 <- ScoreContrib(pc1, obs1 = c(1, 3), obs2 = c(5:10))
plot(Score.Contributions2)
Score.Contributions3 <- ScoreContrib(pc1, obs1 = 1:10, obs2 = 11:15)
plot(Score.Contributions3)
Score.Contributions4 <- ScoreContrib(pc1, obs1 = 1, obs2 = 3)
plot(Score.Contributions4)
}