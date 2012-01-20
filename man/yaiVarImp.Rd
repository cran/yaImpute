\name{yaiVarImp}
\alias{yaiVarImp}
\title{Reprots or plots importance scores for method randomForest}

\description{
  When method \code{randomforest} is used to build a \code{\link{yai}}
  object, the \code{\link[randomForest]{randomForest}} package computes
  variable importance scores. This function computes a composite of the
  scores and scales them. By default the scores are plotted and scores
  themselves are invisibly returned.
}

\arguments{
  \item{object}{an object of class \code{\link{yai}}}
  \item{nTop}{the \code{nTop} most important variables are plotted (returned); if NA or zero, all are returned}
  \item{plot}{if FALSE, no plotting is done, but the scores are returned.}
  \item{\dots}{passed to boxplot function.}
}

\usage{
yaiVarImp(object, nTop=20, plot=TRUE, \dots)
}

\value{
  A data frame with the rows corresponding to the randomForest
  built for each \emph{Y}-variable and the columns corresponding to the
  \code{nTop} most important Y-variables in sorted order.
}

\examples{
data(MoscowMtStJoe)

# get the basal area by species columns
yba  <- MoscowMtStJoe[,1:17]
ybaB <- whatsMax(yba,nbig=7)  # see help on whatsMax

ba <- cbind(ybaB,TotalBA=MoscowMtStJoe[,18])
x <- MoscowMtStJoe[,37:64]
x <- x[,-(4:5)]
rf <- yai(x=x,y=ba,method="randomForest")

yaiVarImp(rf)

keep=colnames(yaiVarImp(rf,plot=FALSE,nTop=9))

newx <- x[,keep]
rf2 <- yai(x=newx,y=ba,method="randomForest")

yaiVarImp(rf2,col="gray")

compare.yai(rf,rf2)

}

\seealso{
\code{\link{yai}}, \code{\link{yaiRFsummary}}, \code{\link{compare.yai}}
}

\author{
  Nicholas L. Crookston \email{ncrookston.fs@gmail.com} \cr
  Andrew O. Finley \email{finleya@msu.edu}
}

\keyword{misc}
\keyword{multivariate}
\keyword{tree}


