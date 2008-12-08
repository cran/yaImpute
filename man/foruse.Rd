\name{foruse}
\alias{foruse}
\title{Report a complete imputation}

\description{
  Provides a matrix of all observations with the reference observation
  identification best used to represent it, followed by the distance.
}

\usage{
foruse(object,kth=NULL)
}

\arguments{
  \item{object}{an object created by \code{\link{yai}}}
  \item{kth}{when NULL, the best pick is reported, otherwise the kth neighbor is picked.}
}

\examples{
require(yaImpute)

data(iris)

# form some test data
refs=sample(rownames(iris),50)
x <- iris[,1:3]      # Sepal.Length Sepal.Width Petal.Length
y <- iris[refs,4:5]  # Petal.Width Species

# build a yai object using mahalanobis
mal <- yai(x=x,y=y,method="mahalanobis",k=3)

foruse(mal)  # for references, use is equal to the rowname
foruse(mal,kth=1) # for references, use is an row to the kth reference.

# get all the choices:
cbind(foruse(mal),foruse(mal,kth=1),foruse(mal,kth=2),foruse(mal,kth=3))

}

\author{
  Nicholas L. Crookston \email{ncrookston@fs.fed.us} \cr
  Andrew O. Finley \email{finleya@msu.edu}
}

\keyword{multivariate}

