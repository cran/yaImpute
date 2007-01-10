# Creates a summary various aspects of several randomForests (one for each variable)
# when method randomForest is used.


yaiRFsummary = function(object, nTop=0)
{
   if (class(object) != "yai") stop ("arg must be of class yai")
   if (object$method != "randomForest") stop ("method must be randomForest")
   if (!require (randomForest)) stop("install randomForest and try again")
   MeanDecreaseGini = yaiVarImp(object, nTop, plot=FALSE)

   oobError  = vector(mode="numeric",length=length(names(object$ranForest)))
   levels    = vector(mode="integer",length=length(names(object$ranForest)))
   ntree     = vector(mode="integer",length=length(names(object$ranForest)))

   i = 0
   for (Rf in object$ranForest)
   {
     i = i+1
     oobError[i] = Rf$err.rate[Rf$ntree,"OOB"]
     levels  [i] = nrow(Rf$confusion)
     ntree   [i] = Rf$ntree
   }
   forestAttributes=data.frame(ntree,oobError,levels)
   rownames(forestAttributes)=names(object$ranForest)
   list(forestAttributes=forestAttributes,MeanDecreaseGini=MeanDecreaseGini)
}

