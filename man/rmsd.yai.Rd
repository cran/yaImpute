\name{rmsd.yai}
\alias{rmsd.yai}

\title{Root Mean Square Difference between observed and imputed}

\description{Computes the root mean square difference (RMSD) between
  observed and imputed values for each observation that has both. RMSD
  is computationally like RMSE, but they differ in interpretation. By default
  the RMSD values are scaled by the standard deviation of the respective variable.
}

\usage{
rmsd.yai (object,vars=NULL,scale=FALSE,\dots)
}

\arguments{
  \item{object}{an object created by \code{\link{yai}} or \code{\link{impute.yai}}}
  \item{vars}{a list of variable names you want to include, if NULL all available
    variables are included}
  \item{scale}{when \code{TRUE}, the values are scaled (see details)}
  \item{\dots}{passed to called methods (not currently used)}
}

\value{
  A data frame with the row names as vars and the column as \code{rmsd}. When
  \code{scale=TRUE}, the column name is \code{rmsdS}.
}

\details{
  By default, RMSD is computed using standard formula for its related statistic,
  RMSE. When \code{scale=TRUE}, RMSD is divided by the standard deviation of the
  observed data. The standard deviation is computed over the original list of
  \emph{reference} observations and is therefore not influenced
  by the number of times a given reference is used for imputation.
}

\seealso{
\code{\link{yai}}, \code{\link{impute.yai}} and \url{./../doc/yaImputePaper.pdf}
(or this alternate \url{http://forest.moscowfsl.wsu.edu/gems/yaImputePaper.pdf}).
}

\author{
  Nicholas L. Crookston \email{ncrookston@fs.fed.us} \cr
  Andrew O. Finley \email{finleya@msu.edu}
}

\keyword{multivariate}


