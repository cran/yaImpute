addXlevels = function (object,origDataFrame)
{
   if (class(object) != "randomForest") stop ("object must be class randomForest")
   if (class(origDataFrame) != "data.frame") stop ("origDataFrame must be class data.frame")
   xlevels=vector(mode="list",length=nrow(object$importance))
   names(xlevels)=rownames(object$importance)
   if (length(intersect(names(xlevels),names(origDataFrame))) != length(names(xlevels)))
      stop ("missing columns in origDataFrame")
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
