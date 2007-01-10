\name{MoscowMtStJoe}
\docType{data}
\alias{MoscowMtStJoe}

\title{Moscow Mountain and St. Joe Woodlands (Idaho, USA) Tree and LiDAR Data}

\description{
  Data used to compare the utility of discrete-return light detection
  and ranging (LiDAR) data and multispectral satellite imagery, and their
  integration, for modeling and mapping basal area and tree density across
  two diverse coniferous forest landscapes in north-central Idaho, USA.
}

\usage{
data(MoscowMtStJoe)
}

\format{
  A data frame with 165 rows and 64 columns: \cr\cr
  Ground based measurements of trees:
  \describe{
    \item \emph{ABGR\_BA } Basal area \eqn{(m^2/ha)} of ABGR
    \item \emph{ABLA\_BA } Basal area \eqn{(m^2/ha)} of ABLA
    \item \emph{ACGL\_BA } Basal area \eqn{(m^2/ha)} of ACGL
    \item \emph{BEOC\_BA } Basal area \eqn{(m^2/ha)} of BEOC
    \item \emph{LAOC\_BA } Basal area \eqn{(m^2/ha)} of LAOC
    \item \emph{PICO\_BA } Basal area \eqn{(m^2/ha)} of PICO
    \item \emph{PIEN\_BA } Basal area \eqn{(m^2/ha)} of PIEN
    \item \emph{PIMO\_BA } Basal area \eqn{(m^2/ha)} of PIMO
    \item \emph{PIPO\_BA } Basal area \eqn{(m^2/ha)} of PIPO
    \item \emph{POBA\_BA } Basal area \eqn{(m^2/ha)} of POBA
    \item \emph{POTR\_BA } Basal area \eqn{(m^2/ha)} of POTR
    \item \emph{PSME\_BA } Basal area \eqn{(m^2/ha)} of PSME
    \item \emph{SAEX\_BA } Basal area \eqn{(m^2/ha)} of SAEX
    \item \emph{THPL\_BA } Basal area \eqn{(m^2/ha)} of THPL
    \item \emph{TSHE\_BA } Basal area \eqn{(m^2/ha)} of TSHE
    \item \emph{TSME\_BA } Basal area \eqn{(m^2/ha)} of TSME
    \item \emph{UNKN\_BA } Basal area \eqn{(m^2/ha)} of unknown species
    \item \emph{Total\_BA} Basal area \eqn{(m^2/ha)} total over all species
    \item \emph{ABGR\_TD } Trees per ha of ABGR
    \item \emph{ABLA\_TD } Trees per ha of ABLA
    \item \emph{ACGL\_TD } Trees per ha of ACGL
    \item \emph{BEOC\_TD } Trees per ha of BEOC
    \item \emph{LAOC\_TD } Trees per ha of LAOC
    \item \emph{PICO\_TD } Trees per ha of PICO
    \item \emph{PIEN\_TD } Trees per ha of PIEN
    \item \emph{PIMO\_TD } Trees per ha of PIMO
    \item \emph{PIPO\_TD } Trees per ha of PIPO
    \item \emph{POBA\_TD } Trees per ha of POBA
    \item \emph{POTR\_TD } Trees per ha of POTR
    \item \emph{PSME\_TD } Trees per ha of PSME
    \item \emph{SAEX\_TD } Trees per ha of SAEX
    \item \emph{THPL\_TD } Trees per ha of THPL
    \item \emph{TSHE\_TD } Trees per ha of TSHE
    \item \emph{TSME\_TD } Trees per ha of TSME
    \item \emph{UNKN\_TD } Trees per ha of unknown species
    \item \emph{Total\_TD} Trees per ha total over all species
  }
  Geographic Location, Slope and Aspect:
  \describe{
    \item \emph{EASTING  } UTM (Zone 11) easting at plot center
    \item \emph{NORTHING } UTM (Zone 11) northing at plot center
    \item \emph{ELEVATION} Mean elevation (m) above sea level over plot
    \item \emph{SLPMEAN  } Mean slope (percent) over plot
    \item \emph{ASPMEAN  } Mean aspect (degrees) over plot
  }
  Advanced Land Imager (ALI):
  \describe{
    \item \emph{B1MEAN } Mean of 30 m ALI band 1 pixels intersecting plot
    \item \emph{B2MEAN } Mean of 30 m ALI band 2 pixels intersecting plot
    \item \emph{B3MEAN } Mean of 30 m ALI band 3 pixels intersecting plot
    \item \emph{B4MEAN } Mean of 30 m ALI band 4 pixels intersecting plot
    \item \emph{B5MEAN } Mean of 30 m ALI band 5 pixels intersecting plot
    \item \emph{B6MEAN } Mean of 30 m ALI band 6 pixels intersecting plot
    \item \emph{B7MEAN } Mean of 30 m ALI band 7 pixels intersecting plot
    \item \emph{B8MEAN } Mean of 30 m ALI band 8 pixels intersecting plot
    \item \emph{B9MEAN } Mean of 30 m ALI band 9 pixels intersecting plot
    \item \emph{PANMEAN} Mean of 10 m PAN band pixels intersecting plot
    \item \emph{PANSTD } Standard deviation of 10 m PAN band pixels intersecting plot
  }
  LiDAR Intensity:
 \describe{
    \item \emph{INTMEAN} Mean of 2 m intensity pixels intersecting plot
    \item \emph{INTSTD } Standard deviation of 2 m intensity pixels intersecting plot
    \item \emph{INTMIN } Minimum of 2 m intensity pixels intersecting plot
    \item \emph{INTMAX } Maximum of 2 m intensity pixels intersecting plot
  }
  LiDAR Height:
  \describe{
    \item \emph{HTMEAN} Mean of 6 m height pixels intersecting plot
    \item \emph{HTSTD } Standard deviation of 6 m height pixels intersecting plot
    \item \emph{HTMIN } Minimum of 6 m height pixels intersecting plot
    \item \emph{HTMAX } Maximum of 6 m height pixels intersecting plot
  }
  LiDAR Canopy Cover:
  \describe{
    \item \emph{CCMEAN} Mean of 6 m canopy cover pixels intersecting plot
    \item \emph{CCSTD } Standard deviation of 6 m canopy cover pixels intersecting plot
    \item \emph{CCMIN } Minimum of 6 m canopy cover pixels intersecting plot
    \item \emph{CCMAX } Maximum of 6 m canopy cover pixels intersecting plot
  }
}

\references{
  Hudak, A.T.; Crookston, N.L.; Evans, J.S.; Falkowski, M.J.; Smith,
  A.M.S.; Gessler, P.E.; Morgan, P. (2006). Regression modeling and mapping
  of coniferous forest basal area and tree density from discrete-return
  lidar and multispectral satellite data. \emph{Can. J. Remote Sensing}.
  32(2):126-138. \url{http://www.treesearch.fs.fed.us/pubs/24612}
}

\source{
  Dr. Andrew T. Hudak  \cr
  USDA Forest Service \cr
  Rocky Mountain Research Station \cr
  1221 South Main \cr
  Moscow, Idaho, USA 83843
}


\keyword{datasets}

