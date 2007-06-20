\name{AsciiGridImpute}
\alias{AsciiGridImpute}
\alias{AsciiGridPredict}
\title{Imputes/Predicts data for Ascii Grid maps}

\description{
   \code{AsciiGridImpute} finds nearest neighbor \emph{reference}
   observations for each point in the input grid maps and outputs maps
   of selected Y-variables in a set of output grid maps. 

   \code{AsciiGridPredict} applies a predict function to each point in the
   input grid maps and outputs maps of the prediction(s) in one or more
   output grid maps (see Details).

  One row of the each grid maps is read and processed at a time thereby
  avoiding the need to build huge objects in R that would be necessary all
  the rows of all the maps were processed together.}

\usage{
AsciiGridImpute(object,xfiles,outfiles,xtypes=NULL,ancillaryData=NULL,
                ann=NULL,lon=NULL,lat=NULL,rows=NULL,cols=NULL,
                nodata=NULL,myPredFunc=NULL,\dots)

AsciiGridPredict(object,xfiles,outfiles,xtypes=NULL,lon=NULL,lat=NULL,
                 rows=NULL,cols=NULL,nodata=NULL,myPredFunc=NULL,\dots)
}

\arguments{
  \item{object}{An object of class \code{\link{yai}}, any other object for which
       a \code{\link[stats]{predict}} function is defined. See Details.}
  \item{xfiles}{A \code{\link[base]{list}} of input file names where there is one
       grid file for each X-variable. List elements must be given the same names
       as the X-variables they correspond with and there must be one file for
       each X-variable used when \code{object} was built.}
  \item{outfiles}{One of these two forms:
      \describe{
       \item (1) A file name that is understood to correspond to the single prediction
       returned by \code{myPredFunc} (generally applies to \code{AsciiGridPredict}), or
       \item (2) A \code{\link[base]{list}} of output file names where there is one
       grid file for each \emph{desired} output variable. While there may be many Y-variables
       in predicted for \code{object}, only those for which an output grid is desire need to be specified.
       In addition to variables names, the following two special names can be coded when the
       object class is \code{yai}: For \code{distance="}\emph{filename}\code{"} a map of the
       distances is output and if \code{useid="}\emph{filename}\code{"} a map of useid's
       is output. When \code{myPredFunc} is not \code{NULL} and when it returns a vector,
       an additional special name of \code{predict="}\emph{filename}\code{"} is used.} }
  \item{xtypes}{A list of data types names that corresponds exactly to data type of the
       maps listed in \code{xfiles}. Each value can be one of:
       \code{"logical", "integer", "numeric", "character"}. If NULL,
       or if a type is missing for a member of \code{xfiles}, type \code{"numeric"} is used.}
  \item{ancillaryData}{A data frame of Y-variables that may not have been used in
       the original call to \code{\link{yai}}. There must be one row for
       each reference observation, no missing data, and row names must match those used
       in the original reference observations.}
  \item{ann}{if NULL, the value is taken from \code{object}. When TRUE, \code{\link{ann}} is
       used to find neighbors, and when FALSE a slow exact search is used (ignored for when 
       method randomForest is used when the original \code{\link{yai}} object was created).}
  \item{lon}{if NULL, the value of \code{rows} is used. Otherwise, a 2-element
       vector given the range of longitudes (vertical distance) desired for the output.} 
  \item{lat}{if NULL, the value of \code{cols} is used. Otherwise, a 2-element
       vector given the range of latitudes (horizontal distance) desired for the output.} 
  \item{rows}{if NULL, all rows from the input grids are used. Otherwise, rows is a 2-element
       vector given the rows desired for the output. If the second element is greater than
       the number of rows, the header value \code{YLLCORNER} is adjusted accordingly. Ignored
       if \code{lon} is specified.}
  \item{cols}{if NULL, all columns from the input grids are used. Otherwise, cols is a 2-element
       vector given the columns desired for the output. If the first element is greater than
       one, the header value \code{XLLCORNER} is adjusted accordingly. Ignored
       if \code{lat} is specified.}
  \item{nodata}{the \code{NODATA_VALUE} for the output. If NULL, the value is taken from the
       input grids.}
  \item{myPredFunc}{called to predict output using the \code{object} and newdata from
       the \code{xfiles}. Two arguments are required, the first is \code{object} and the
       second is a data frame of the new predictor variables. If NULL, the generic
       \code{predict} function is called for \code{object}.}
  \item{\dots}{passed to \code{myPredFunc} or \code{predict}.}
}

\details{
  The input maps are assumed to be Ascii Grid maps with 6-line headers
  containing the following tags: \code{NCOLS, NROWS, XLLCORNER, YLLCORNER,
  CELLSIZE} and \code{NODATA_VALUE} (case insensitive).

  The function builds data frames from the input maps (one for each row) then 
  builds predictions using those data frames as new data. Each row of the 
  input maps is processed in sequence so that the entire maps are not stored 
  in memory. The function works by opening all the input and reads one line 
  (row) at a time from each. The output file(s) are created one line at time 
  as the input maps are processed. Missing values are predicted when a 
  missing value is in the input (that is, the predict function is not called 
  for any grid cell that contains a missing value and the missing value set
  for the output grid cell automatically). 
  
  Use \code{AsciiGridImpute} for objects builts with \code{\link{yai}}, 
  otherwise use \code{AsciiGridPredict}. When \code{AsciiGridPredict} is 
  used, the following rules apply. First, when \code{myPredFunc} is not 
  null it is called with the arguments \code{object, newdata} where the 
  new data is the data frame built from the input maps, otherwise the 
  generic \code{\link{predict}} function is called with these same arguments. 
  When \code{object} and \code{myPredFunc} are both NULL a copy 
  \code{newdata} used as the prediction. This is useful when \code{lat, lon, rows,} 
  or \code{cols} are used in to subset the maps.

  If factors are used as X-variables in \code{object}, the levels found the 
  new data are checked against those used in building the \code{object}. If 
  new levels are found, the corresponding output map grid point is set to a 
  missing value (most predict functions will fail in these 
  circumstances). Checking on factors depends on \code{object} containing a 
  meaningful member named \code{xlevels}, as done objects model objects 
  produced by \code{\link[stats]{lm}}. Note that objects produced by 
  \code{\link[randomForest]{randomForest}} version 4.5-19 and prior do 
  not contain \code{xlevels}; use function \code{\link{addXlevels}} to 
  add the necessary element if one or more 
  \code{\link[base]{factor}}\code{s} are used. 
}  
  

\seealso{
\code{\link{yai}}, \code{\link{impute}}, and \code{\link{newtargets}}
}

\value{
A list of factor levels that are in the maps, but not used to build the \code{object}, or
NULL if there are none. The list includes the number of map grid points the offending
level was found.
}

\examples{

## These commands write new files to your working directory

# Use the iris data
data(iris)

# Change species from a character factor to numeric (the sp classes
# can not handle character data).

iris[,5]<-as.factor(as.numeric(iris[,5]))
sLen <- matrix(iris[,1],10,15)
sWid <- matrix(iris[,2],10,15)
pLen <- matrix(iris[,3],10,15)
pWid <- matrix(iris[,4],10,15)
spcd <- matrix(as.numeric(iris[,5]),10,15)

# Make a "map" of each variable.

header = c("NCOLS 15","NROWS 10","XLLCORNER 1","YLLCORNER 1",
           "CELLSIZE 1","NODATA_VALUE -9999")
cat(file="slen.txt",header,sep="\n")
cat(file="swid.txt",header,sep="\n")
cat(file="plen.txt",header,sep="\n")
cat(file="pwid.txt",header,sep="\n")
cat(file="spcd.txt",header,sep="\n")


write.table(sLen,file="slen.txt",append=TRUE,col.names=FALSE,
            row.names=FALSE)
write.table(sWid,file="swid.txt",append=TRUE,col.names=FALSE,
            row.names=FALSE)
write.table(pLen,file="plen.txt",append=TRUE,col.names=FALSE,
            row.names=FALSE)
write.table(pWid,file="pwid.txt",append=TRUE,col.names=FALSE,
            row.names=FALSE)
write.table(spcd,file="spcd.txt",append=TRUE,col.names=FALSE,
            row.names=FALSE)

# sample the data
refs <- sample(rownames(iris),50)
y <- data.frame(Species=iris[refs,5],row.names=rownames(iris[refs,]))

# build a yai imputation for the reference data.
rfNN <- yai(x=iris[refs,1:4],y=y,method="randomForest")

xfiles <- list(Sepal.Length="slen.txt",Sepal.Width="swid.txt",
               Petal.Length="plen.txt",Petal.Width="pwid.txt")
outfiles1 <- list(distance="dist.txt",Species="spOutrfNN.txt")

# map the imputation-based predictions for the input maps
AsciiGridImpute(rfNN,xfiles,outfiles1,ancillaryData=iris)

# build a randomForest predictor
rf <- randomForest(x=iris[refs,1:4],y=iris[refs,5])

# map the predictions for the input maps
outfiles2 <- list(predict="spOutrf.txt")
AsciiGridPredict(rf,xfiles,outfiles2,xtypes=NULL,rows=NULL)

# read the species maps and plot them using the sp package classes
if (require(sp)) {
   spOrig <- read.asciigrid("spcd.txt")
   sprfNN <- read.asciigrid("spOutrfNN.txt")
   sprf <- read.asciigrid("spOutrf.txt")
   dist <- read.asciigrid("dist.txt")

   par(mfcol=c(2,2))
   image(spOrig,col=c(1,2,3))
   title("Original")
   image(sprfNN,col=c(1,2,3))
   title("Using Predict")
   image(sprf,col=c(1,2,3))
   title("Using Impute")
   image(dist)
   title("Neighbor Distances")
}
}


\author{
  Nicholas L. Crookston \email{ncrookston@fs.fed.us} \cr
  Andrew O. Finley \email{finleya@msu.edu}
}

\keyword{spatial}
\keyword{utilities}
