# Creates a summary of the variable importance scores over several randomForests (one
# for each variable) when method randomForest is used. These values are then scaled.

yaiVarImp = function(object, nTop=20, plot=TRUE, ...)
{
   if (class(object) != "yai") stop ("arg must be of class yai")
   if (object$method != "randomForest") stop ("method must be randomForest")
   if (!require (randomForest)) stop("install randomForest and try again")
   scaledImportance = matrix(NA, nrow = length(names(object$ranForest)), ncol=length(xvars(object)))
   colnames(scaledImportance) = xvars(object)
   rownames(scaledImportance) = names(object$ranForest)

   i = 0
   for (Rf in object$ranForest)
   {
     i = i+1
     one = importance(Rf)
     scaledImportance[i,] = scale(one[,ncol(one)])
   }

   if (is.na(nTop) | nTop == 0) nTop=ncol(scaledImportance)
   scaledImportance = data.frame(scaledImportance)
   nTop = min(ncol(scaledImportance), nTop)
   best = sort(apply(scaledImportance, 2, median), decreasing = TRUE, index.return = TRUE)$ix[1:nTop]
   if (plot)
   {
      plt = par()$plt
      oldplt = plt
      plt[1] = .2
      boxplot(as.data.frame(scaledImportance[,best]), horizontal=TRUE, par(plt=plt), las=1,
              main=deparse(substitute(object)), xlab="Scaled Importance",...)
      par(plt=oldplt)
      invisible(scaledImportance[,best,FALSE])
   }
   else return(scaledImportance[,best,FALSE])
}

