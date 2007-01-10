\name{plot.compare.yai}
\alias{plot.compare.yai}
\title{Plots a compare.yai object}

\description{
  Provides a matrix of plots for objects created by \code{\link{compare.yai}}.
}

\usage{
\method{plot}{compare.yai}(x,\dots,colors=1,lineColor=2)
}

\arguments{
  \item{x}{a data frame created by \code{\link{compare.yai}}.}
  \item{\dots}{passed to plot functions.}
  \item{colors}{the colors (recycled as necessary) used for the points.}
  \item{lineColor}{the color of the 1:1 line.}
}

\seealso{
\code{\link{yai}}, \code{\link{compare.yai}},
\code{\link{impute.yai}}, \code{\link{rmsd.yai}}}
}

\author{
  Nicholas L. Crookston \email{ncrookston@fs.fed.us} \cr
  Andrew O. Finley \email{afinley@stat.umn.edu}
}

\keyword{hplot}


