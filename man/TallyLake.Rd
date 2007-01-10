\name{TallyLake}
\docType{data}
\alias{TallyLake}

\title{Tally Lake, Flathead National Forest, Montana, USA}

\description{
  Polygon-based reference data used by Stage and Crookston (In press)
  to demonstrate partitioning of error components and related statistics.
  Observations are summaries of data collected on forest stands (ploygons).
}

\usage{
data(TallyLake)
}

\format{
  A data frame with 847 rows and 29 columns: \cr\cr
  Ground based measurements of trees (Y-variables):
  \describe{
     \item \emph{TopHt     } Height of tallest trees (ft)
     \item \emph{LnVolL    } Log of the volume \eqn{(ft^3/acre)} of western larch
     \item \emph{LnVolDF   } Log of the volume \eqn{(ft^3/acre)} of Douglas\-fir
     \item \emph{LnVolLP   } Log of the volume \eqn{(ft^3/acre)} of lodgepole pine
     \item \emph{LnVolES   } Log of the volume \eqn{(ft^3/acre)} of Engelmann spruce
     \item \emph{LnVolAF   } Log of the volume \eqn{(ft^3/acre)} of alpine fir
     \item \emph{LnVolPP   } Log of the volume \eqn{(ft^3/acre)} of ponderosa pine
     \item \emph{CCover    } Canopy cover (percent)
  }
  Geographic Location, Slope, and Aspect (X-variables):
  \describe{
    \item \emph{utmx      } UTM easting at plot center
    \item \emph{utmy      } UTM northing at plot center
    \item \emph{elevm     } Mean elevation (ft) above sea level over plot
    \item \emph{eevsqrd   } \eqn{(elevm-1600)^2}
    \item \emph{slopem    } Mean slope (percent) over plot
    \item \emph{slpcosaspm} Mean of slope (proportion) times the cosine of aspect (see
     Stage (1976) for description of this transformation)
    \item \emph{slpsinaspm} Mean of slope (proportion) times the sine of aspect
  }
  Additional X-variables:
  \describe{
    \item \emph{ctim    } Mean of slope curviture over pixels in stand
    \item \emph{tmb1m   } Mean of LandSat band 1 over pixels in stand
    \item \emph{tmb2m   } Mean of LandSat band 2 over pixels in stand
    \item \emph{tmb3m   } Mean of LandSat band 3 over pixels in stand
    \item \emph{tmb4m   } Mean of LandSat band 4 over pixels in stand
    \item \emph{tmb5m   } Mean of LandSat band 5 over pixels in stand
    \item \emph{tmb6m   } Mean of LandSat band 6 over pixels in stand
    \item \emph{durm    } Mean of light duration over pixels in stand
    \item \emph{insom   } Mean of solar insolation over pixels in stand
    \item \emph{msavim  } Mean of AVI for pixels in stand
    \item \emph{ndvim   } Mean of NDVI for pixels in stand
    \item \emph{crvm    } Mean of slope curviture for pixels in stand
    \item \emph{tancrvm } Mean of tangent curvature for pixels in stand
    \item \emph{tancrvsd} Standard deviation of tangent curvature for pixels in stand
  }
}

\references{
  Stage, A.R.; Crookston, N.L. (In press). Partitioning error components
  for accuracy-assessment of near neighbor methods of imputation.
  \emph{For. Sci.} \url{./../doc/StagePartitioningFS.pdf}

  Stage, A.R. (1976). An expression for the effect of aspect, slope, and
  habitat type on tree growth. \emph{For. Sci.} 22(4):457-460.
}

\source{
  USDA Forest Service
}


\keyword{datasets}

