# Imputes the observation for a variable from a reference observation to a
# target observation. Also, imputes a value for a reference from other
# references. This practice is useful for validation.
#
# Arguments:
#   object is of class yai (built using function yai)
#   ancillaryData is a data frame of reference, target, or both observations
#       with variables that may or may not be in original problem.
#       Imputations are made for these variables using the neighbor
#       relationships found in object.
#   method defines how continuous variables are imputed, where:
#       closest     = the value from the nearest neighbor is imputed for
#                     target observations. For reference observations,
#                     the value form the observation itself is excluded.
#                     When k==1, closest is always used. 
#       allk        = a column for each variable and for each k is returned. 
#       onek        = a column for each variable and for the kth neighbor is 
#                     returned.
#       mean        = An average over the k neighbors is taken
#       dstWeighted = a weighted average is taken over the k
#                     neighbors where the weights are 1/(1+d)
#  method.factor defines how factors are imputed, where: 
#       if method is allk or onek, method.factor is ignored.
#       closest     = same a continous (always used with k==1).
#       random      = random pick with sampling proportional to 1/max(d,.0001)
#                     NOT YET SUPPORTED
#       overallfactors = the closest considering all factors
#                     NOT YET SUPPORTED
#  k is the number of neighbors to use. When NULL, all those available
#    in the yai object are used.
#  vars a character vector of variable names desired. If NULL (the default),
#    then the list is formed by taking all variables in the input object.
#  observed is a flag, when TRUE, the returned data frame includes
#    the observed values (see value, below)
# Value:
#  A data frame of class impute.yai with all variables
#  in the problem or just those listed in "vars" if it is not
#  null. When observed is TRUE, the set of variables
#  is duplicated with .o added to the names and each is
#  the observed value when one exists for the observation.

impute <- function(object,...) UseMethod("impute")


impute.yai <- function (object,ancillaryData=NULL,method="closest",
                       method.factor="closest",k=NULL,
                       vars=NULL,observed=TRUE,...)
{
   # ======== predict functions used internally
   findFactors =  get("findFactors",asNamespace("yaImpute"))

   pred.c <- function (refs,ids,w=NULL,method="closest",k=1,vars,observed)
   {
      if (is.null(vars)) vars <- colnames(refs)
      else vars <- intersect(vars,colnames(refs))
      if (is.null(vars) | length(vars)==0) return (NULL)

      if (method=="closest" | k==1)
      {
         if (length(vars)==1)
         {
            ans <- as.matrix(refs[ids[,1],vars])
            colnames(ans) <- vars
         }
         else ans <- refs[ids[,1],vars]
      }
      else
      {
         ans <- matrix(data=NA, nrow = nrow(ids), ncol = length(vars))
         colnames(ans) <- vars
      }
      rownames(ans) <- rownames(ids)
      ans <- as.data.frame(ans)

      if (method=="mean")
      {
         for (row in 1:nrow(ans)) ans[row,vars] <- apply(refs[ids[row,1:k],vars],2,mean)
      }
      if (method=="dstWeighted")
      {
         if (is.null(w))
         {
            warning ("w is required when method is dstWeighted")
            return(NULL)
         }
         for (row in 1:nrow(ans))
         {
            wei <- 1/(1+w[row,1:k])
            wei/sum(wei)
            ans[row,vars] <- apply(refs[ids[row,1:k],vars],2,weighted.mean,wei)
         }
      }
      if (observed)
      {
         obs <- matrix(data=NA, nrow = nrow(ans), ncol = ncol(ans))
         rownames(obs) <- rownames(ans)
         colnames(obs) <- vars
         obs <- as.data.frame(obs)
         commonRows <- intersect(rownames(ans),rownames(refs))
         if (length(commonRows)>0) obs[commonRows,vars] <- refs[commonRows,vars]
         colnames(obs) <- paste(vars,"o",sep=".")
         ans <- cbind(ans,obs)
      }
      ans
   }
   pred.f <- function (refs,ids,w=NULL,method="closest",k=1,vars,observed)
   {
      if (is.null(vars)) vars <- colnames(refs)
      else vars <- intersect(vars,colnames(refs))
      if (is.null(vars) | length(vars)==0) return (NULL)
      if (method != "closest")
      {
         warning ("only method closest is supported for factors")
         method <- "closest"
      }
      if (method=="closest" | k==1)
      {
         if (length(vars)==1)
         {
            ans <- data.frame(refs[ids[,1],vars])
            colnames(ans) <- vars
         }
         else ans <- data.frame(refs[ids[,1],vars])
      }
      rownames(ans) <- rownames(ids)
      if (observed)
      {
         obs <- matrix(data=NA, nrow = nrow(ans), ncol = ncol(ans))
         rownames(obs) <- rownames(ans)
         colnames(obs) <- vars
         obs <- as.data.frame(obs)
         commonRows <- intersect(rownames(ans),rownames(refs))
         if (length(commonRows)>0)
         {
            for (var in vars)
            {
               obs[commonRows,var] <- levels(refs[,var])[refs[commonRows,var]]
               obs[,var] <- factor(obs[,var])
           }
         }
         names(obs) <- paste(vars,"o",sep=".")
         ans <- cbind(ans,obs)
      }
      ans
   }
   pred <- function (refs,ids,w=NULL,method="closest",
                    method.factor="closest",k=1,vars,observed)
   {
      factors <- findFactors(refs)
      nfactors <- sum(factors)
      if (nfactors>0 && method.factor != "closest" && k==1)
      {
         warning ("method.factor was set to closest because k=1")
         method.factor <- "closest"
      }
      if (nfactors == 0)
         out <- pred.c(refs=refs,ids=ids,w=w,method=method,
                    k=k,vars=vars,observed=observed)
      else if (nfactors == ncol(refs))
         out <- pred.f(refs=refs,ids=ids,w=w,method=method.factor,
                    k=k,vars=vars,observed=observed)
      else
      {
         tmp <- data.frame(refs[,!factors],row.names=rownames(refs))
         colnames(tmp) <- colnames(refs)[!factors]
         p1 <- pred.c(refs=tmp,ids=ids,w=w,method=method,
                   k=k,vars=vars,observed=observed)
         tmp <- data.frame(refs[,factors],row.names=rownames(refs))
         colnames(tmp) <- colnames(refs)[factors]
         p2 <- pred.f(refs=tmp,ids=ids,w=w,method=method.factor,
                   k=k,vars=vars,observed=observed)
         if      (is.null(p1) && is.null(p2)) out <- NULL
         else if (is.null(p1)) out <- p2
         else if (is.null(p2)) out <- p1
         else                  out <- cbind(p1,p2)
      }
      out
   }

   # ===========================


   if (missing(object)) stop ("object required.")
   if (class(object) != "yai") stop ("class must be yai")

   if (is.null(vars))
   {
      if (is.null(ancillaryData)) 
      {
         if (object$method != "randomForest") vars <- yvars(object)
         else if (names(object$ranForest)[[1]] == "unsupervised") vars <- xvars(object)
      }            
      else vars <- colnames(ancillaryData)
   }
   posMethods <- c("closest","mean","dstWeighted")
   if (length(intersect(method,posMethods))==0)
      stop (paste("method=",method," must be one of: {",
            paste(posMethods,collapse=", "),"}",sep=""))

   if (is.null(k)) k <- object$k
   if (k>object$k | k==0)
   {
      warning ("k out of range, set to ",object$k)
      k <- object$k
   }

   if (method != "closest" && k==1)
   {
      warning ("method was set to closest because k==1")
      method <- "closest"
   }

   if (is.null(ancillaryData))
   {
      r <- NULL
      if (length(object$neiIdsRefs)>0)
      {
         if (!(ncol(object$yRefs) == 1 && names(object$yRefs)[1]=="ydummy"))
              yPredRefs <- pred(refs=object$yRefs,ids=object$neiIdsRefs,w=object$neiDstRefs,
                           method=method,method.factor=method.factor,k=k,vars=vars,observed=observed)
         else yPredRefs <- NULL
         xPredRefs <- pred(refs=object$xall,ids=object$neiIdsRefs,w=object$neiDstRefs,
                           method=method,method.factor=method.factor,k=k,vars=vars,observed=observed)
         if      (is.null(yPredRefs) && is.null(xPredRefs)) r <- NULL
         else if (is.null(yPredRefs)) r <- xPredRefs
         else if (is.null(xPredRefs)) r <- yPredRefs
         else                         r <- cbind(yPredRefs,xPredRefs)

      }
      t <- NULL
      if (length(object$neiIdsTrgs)>0)
      {
         if (!(ncol(object$yRefs) == 1 && names(object$yRefs)[1]=="ydummy"))
             yPredTrgs <- pred(refs=object$yRefs,ids=object$neiIdsTrgs,w=object$neiDstTrgs,
                               method=method,method.factor=method.factor,k=k,vars=vars,observed=observed)
         else yPredTrgs <- NULL
         xPredTrgs <- pred(refs=object$xall,ids=object$neiIdsTrgs,w=object$neiDstTrgs,
                           method=method,method.factor=method.factor,k=k,vars=vars,observed=observed)

         if      (is.null(yPredTrgs) && is.null(xPredTrgs)) t <- NULL
         else if (is.null(yPredTrgs)) t <- xPredTrgs
         else if (is.null(xPredTrgs)) t <- yPredTrgs
         else                         t <- cbind(yPredTrgs,xPredTrgs)
      }
      if      (is.null(r) && is.null(t)) out <- NULL
      else if (is.null(r)) out <- t
      else if (is.null(t)) out <- r
      else                 out <- rbind(r,t)

      scale <- data.frame(center=c(object$xScale$center,object$yScale$center),
                           scale=c(object$xScale$scale, object$yScale$scale))
   }
   else
   {
      out <- NULL
      if (!is.null(vars))
      {
         ancillaryData=ancillaryData[,vars,FALSE]
         if (is.null(ncol(ancillaryData))) stop ("requested variables not present in ancillaryData")
      }
      rownames(ancillaryData)=as.character(rownames(ancillaryData))
      ids <- as.character(rownames(object$xRefs))
      common <- intersect(ids,rownames(ancillaryData))
      missing <- setdiff(common,rownames(ancillaryData))
      if (length(missing) != 0) warning (paste("no data for",length(missing),"observations:",
                                paste(missing[1:min(15,length(missing))],collapse=",")))
      w  <- rbind(object$neiDstRefs,object$neiDstTrgs)

      ids <- rbind(object$neiIdsRefs,object$neiIdsTrgs)
      out <- pred(refs=ancillaryData,ids=ids,w=w,method=method,
                  method.factor=method.factor,k=k,vars=vars,observed=observed)

      # find the sd/mean of the data and attach these as an attribute

      notFactors <- !findFactors(ancillaryData)
      if (sum(notFactors) > 0)
      {
         scale <- matrix(data=NA,nrow=ncol(ancillaryData),ncol=2)
         rownames(scale) <- colnames(ancillaryData)
         colnames(scale) <- c("center","scale")
         scale[notFactors,"center"] <- mean(ancillaryData[,rownames(scale)[notFactors]],na.rm=TRUE)
         scale[notFactors,"scale" ] <- sd  (ancillaryData[,rownames(scale)[notFactors]],na.rm=TRUE)
      }
      else scale=NULL
   }
   if (!is.null(out))
   {
      class(out) <- c("impute.yai","data.frame")
      if (!is.null(scale)) attr(out,"scale") <- scale
   }
   out
}
