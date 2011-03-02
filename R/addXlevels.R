addXlevels = function (object,origDataFrame)
{
	 .Deprecated(" "," ","Newer versions of package randomForest do not require addXlevels") 
   if (class(object) != "randomForest") stop ("object must be class randomForest")
   if (class(origDataFrame) != "data.frame") stop ("origDataFrame must be class data.frame")
   xlevels=vector("list",length(object$forest$ncat)) 
   names(xlevels)=names(object$forest$ncat)
   if (length(intersect(names(xlevels),names(origDataFrame))) != length(names(xlevels)))
      stop ("Variables used in training data missing from origDataFrame: ",paste(setdiff(names(origDataFrame),names(xlevels)),collapse=","))
   for (var in names(xlevels))
   {
      ivar = match(var,names(origDataFrame))
      if (is.na(ivar)) stop ("a variable used in object is not in origDataFrame")
      if (is.factor(origDataFrame[,ivar])) xlevels[[var]]=levels(origDataFrame[,ivar])
      else xlevels[[var]]=NULL
   }
   if (length(xlevels)>0) object$xlevels=xlevels
   object
}
