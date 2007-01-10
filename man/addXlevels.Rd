\name{addXlevels}
\alias{addXlevels}
\title{Adds xlevels to randomForest objects}

\description{
  This function adds \code{xlevels} (see \code{\link[stats]{lm}}) to
  \code{\link[randomForest]{randomForest}} objects. Function
  \code{\link{AsciiGridImpute}} will check the levels on the input
  maps to insure that only those used in the fitting are used in the predictions.

 }

\usage{
addXlevels(object,origDataFrame)
}

\seealso{
  \code{\link{yai}} and \code{\link{AsciiGridPredict}}
}

\arguments{
  \item{object}{an object built by \code{\link[randomForest]{randomForest}}}
  \item{origDataFrame}{the data frame used in the \code{\link[randomForest]{randomForest}} fit.}
}

\value{
  An object of class \code{randomForest} with \code{xlevels} appended.
}

\examples{

if (require(randomForest))
{
  data(iris)
  rf = randomForest(x=iris[,2:5],y=iris[,1])
  new = addXlevels(rf,iris)
  print(new$xlevels)
  predict(new)
}
}

\author{
  Nicholas L. Crookston \email{ncrookston@fs.fed.us} \cr
  Andrew O. Finley \email{afinley@stat.umn.edu}
}

\keyword{misc}
\keyword{multivariate}


