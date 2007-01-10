\name{plot.yai}
\alias{plot.yai}
\alias{plot.impute.yai}
\title{Plot observed verses imputed data}

\description{
  Provides a matrix of plots of observed verses imputed values for variables in an
  object created by \code{\link{impute.yai}}, which are of class
  \code{c("impute.yai","data.frame")}.
}

\usage{
\method{plot}{yai}(x,\dots,vars=NULL,colxy=NULL,colspine=NULL,residual=FALSE)
}

\arguments{
  \item{x}{(1) a data frame created by \code{\link{impute.yai}}, or (2) an object
    created by \code{\link{yai}}.}
  \item{\dots}{passed to called functions.}
  \item{vars}{a list of variable names you want to include, if NULL all available
    Y-variables are included.}
  \item{colxy}{a color vector for the xy plots (continuous variables).}
  \item{colspine}{a color vector for the spine plots (factors).}
  \item{residual}{plots in a residual format (observed-imputed over imputed).}
}

\examples{
require(yaImpute)

data(iris)

# form some test data
refs=sample(rownames(iris),50)
x <- iris[,1:3]      # Sepal.Length Sepal.Width Petal.Length
y <- iris[refs,4:5]  # Petal.Width Species

mal <- yai(x=x,y=y,method="mahalanobis")
malImp=impute(mal,newdata=iris)
plot(malImp)

}

\author{
  Nicholas L. Crookston \email{ncrookston@fs.fed.us} \cr
  Andrew O. Finley \email{afinley@stat.umn.edu}
}

\keyword{hplot}


