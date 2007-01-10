# Creates a summary of the variable importance scores over several randomForests (one
# for each variable) when method randomForest is used.

yaiVarImp = function(object, nTop=20, plot=TRUE)
{
   if (class(object) != "yai") stop ("arg must be of class yai")
   if (object$method != "randomForest") stop ("method must be randomForest")
   if (!require (randomForest)) stop("install randomForest and try again")
   MeanDecreaseGini = matrix(0, nrow = length(names(object$ranForest)), ncol=length(xvars(object)))
   colnames(MeanDecreaseGini) = xvars(object)
   rownames(MeanDecreaseGini) = names(object$ranForest)

   i = 0
   for (Rf in object$ranForest)
   {
     i = i+1
     MeanDecreaseGini[i,] = importance(Rf)[,"MeanDecreaseGini"]
   }

   if (is.na(nTop) | nTop == 0) nTop=ncol(MeanDecreaseGini)
   MeanDecreaseGini = data.frame(MeanDecreaseGini)
   nTop = min(ncol(MeanDecreaseGini), nTop)
   best = sort(apply(MeanDecreaseGini, 2, median), decreasing = TRUE, index.return = TRUE)$ix[1:nTop]
   if (plot)
   {
      plt = par()$plt
      oldplt = plt
      plt[1] = .2
      boxplot(as.data.frame(MeanDecreaseGini[,best]), horizontal=TRUE, par(plt=plt), las=1,
              main=deparse(substitute(object)), xlab="MeanDecreaseGini")
      par(plt=oldplt)
      invisible(MeanDecreaseGini[,best])
   }
   else return(MeanDecreaseGini[,best])
}

