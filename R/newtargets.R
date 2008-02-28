# Arguments:
#   object is of class "yai" created by function yai.
#   newdata is a matrix or dataframe that contains data for x variables
#   k is the number of neighbors desired (see function yai).
#
# Value:  A list of class "yai", that is a copy of the input object with
#         the following members replaced:
#
#   call: the call
#
#   obsDropped: a list of the rownames for observations dropped for various
#               reasons (missing data).
#
#   trgRows: a list of the rownames for target observations as a subset
#            of observations in xall (normally, all rows).
#
#   xall: the x variables for all observations.
#
#   neiDstTrgs: A data frame of distances between a target (identified by it's
#      rowname) and the k references. There are k columns.
#
#   neiIdsTrgs: A data frame of reference identifications that correspond to
#       neiDstTrgs.
#
#   ann: use ann or not...if null, the value is taken from object.
#

newtargets=function(object,newdata,ann=NULL)
{

   if (class(object) != "yai") stop ("object must be class yai")
   if (is.null(newdata) | nrow(newdata)==0) stop ("newdata is required")
   if (object$method == "gnn") # (GNN), make sure we have package vegan loaded
      if (!require (vegan)) stop("install vegan and try again")
   if (object$method == "randomForest") # make sure we have package randomForest loaded
      if (!require (randomForest)) stop("install randomForest and try again")

   sumSqDiff=function(x,y) { d=x-y; sum(d*d) }
   factorMatch = get("factorMatch",asNamespace("yaImpute"))

   if (is.null(ann)) ann=object$ann

   object$call=match.call()

   obsDropped=NULL

   # don't redo the factor matching for objects that come already done.
   if (!is.null(attr(newdata,"illegalLevelCounts")) &&
        length(intersect(xvars(object),names(object$xlevels))) > 0)
   {
      newdata = factorMatch(newdata,object$xlevels)
      if (is.list(attr(newdata,"illegalLevelCounts")))
      {
         warning ("NA's generated due to illegal level(s).")
         cat ("Illegal levels\n")
         print(attr(newdata,"illegalLevelCounts"))
      }
   }

   if (is.null(object$theFormula))
   {
      vold =colnames(object$xall)
      have=intersect(vold,colnames(newdata))
      if (length(have) != length(vold))
          stop(paste("required column(s) missing:",paste(is.null(have) ? vold : vold[-have],collapse=", ")))
      xall=na.omit(newdata[,have])
      obsDropped=names(attributes(na.omit(xall))$na.action)
      if (length(obsDropped)>0) warning (nrow(newdata)-nrow(xall)," observation(s) removed")
   }
   else
   {
      xall=model.frame(object$theFormula$x,newdata)
      if (!is.null(object$xDrop)) xall=xall[,!object$xDrop,drop=FALSE]
      obsDropped=setdiff(rownames(newdata),rownames(xall))
      if (length(obsDropped)) warning (length(obsDropped)," observation(s) removed")
   }
   if (nrow(xall) == 0) stop ("no observations")
   trgs=setdiff(rownames(xall),rownames(object$xRefs))
   if (nrow(xall) != length(trgs))
   {
      obsDropped=union(obsDropped,intersect(rownames(object$xRefs),rownames(xall)))
      warning (nrow(xall)-length(trgs)," row(s) in newdata are original references and ignored")
   }

   theCols = colnames(object$xRefs)  # may be changed for reduced rank, depending on method.

   if (object$method == "msn" | object$method == "msn2" |
       object$method == "mahalanobis" | object$method == "ica")
   {
      theCols = rownames(object$projector)
      xcvRefs=scale(object$xRefs,center=object$xScale$center,scale=object$xScale$scale)
      if (length(theCols)<ncol(xcvRefs)) xcvRefs=xcvRefs[,theCols,drop=FALSE]
   }

   if (object$method == "euclidean")
      xcvRefs=scale(object$xRefs,center=object$xScale$center,scale=object$xScale$scale)

   xTrgs=as.data.frame(xall[trgs,theCols,drop=FALSE]) # this is needed by gnn unscalled.
   if (nrow(xTrgs)==0) stop("no observations")

   if (object$method == "gnn") # GNN
   {
      # create a projected space for the reference observations
      predCCA = get("modified.predict.cca",asNamespace("yaImpute"))
      xcvRefs=predCCA(object$ccaVegan,type="lc",rank="full")
      xcvRefs=xcvRefs %*% diag(sqrt(object$ccaVegan$CCA$eig/sum(object$ccaVegan$CCA$eig)))

      # create a projected space for the unknowns (target observations)
      xcvTrgs=scale(xTrgs,center=object$xScale$center,scale=object$xScale$scale)
      xcvTrgs=predCCA(object$ccaVegan,newdata=as.data.frame(xcvTrgs),type="lc",rank="full")
      xcvTrgs=xcvTrgs %*% diag(sqrt(object$ccaVegan$CCA$eig/sum(object$ccaVegan$CCA$eig)))
      nVec = ncol(xcvRefs)
   }
   if (object$method == "randomForest") # randomForest
   {
      nodes=NULL
      for (i in 1:length(object$ranForest))
      {
         predRF = getS3method("predict","randomForest")
         nodeset=attr(predRF(object$ranForest[[i]],rbind(object$xRefs,xTrgs),
                      proximity=FALSE,nodes=TRUE),"nodes")
         if (is.null(nodeset)) stop("randomForest did not return nodes")
         colnames(nodeset)=paste(colnames(nodeset),i,sep=".")
         nodes=if (is.null(nodes)) nodeset else cbind(nodes,nodeset)
      }
      INTrefNodes=as.integer(nodes[rownames(object$xRefs),])
      INTnrow=as.integer(nrow(object$xRefs))
      INTncol=as.integer(ncol(nodes))
   }
   else
   {
      if (object$method != "raw") xcvTrgs=scale(xTrgs,center=object$xScale$center,scale=object$xScale$scale)[,theCols]
      else                        xcvTrgs=xTrgs[,theCols]
      if (!is.null(object$projector))
      {
         xcvRefs=xcvRefs[,theCols,drop=FALSE] %*% object$projector
         xcvTrgs=xcvTrgs[,theCols,drop=FALSE] %*% object$projector
      }
   }

   neiDstTrgs=matrix(data=NA,nrow=nrow(xTrgs),ncol=object$k)
   rownames(neiDstTrgs)=rownames(xTrgs)
   colnames(neiDstTrgs)=paste("Dst.k",1:object$k,sep="")
   neiIdsTrgs=neiDstTrgs
   colnames(neiIdsTrgs)=paste("Id.k",1:object$k,sep="")

   if (object$method != "randomForest")
   {
      if (ann & nrow(xcvTrgs)>0)
      {
          k=object$k
          ann.out=ann(xcvRefs, xcvTrgs, k, verbose=FALSE)$knnIndexDist
          neiDstTrgs[TRUE]=sqrt(ann.out[,(k+1):ncol(ann.out)])
          for (i in 1:k)
             neiIdsTrgs[,i]=rownames(xcvRefs)[ann.out[,i]]
          rownames(neiDstTrgs)=rownames(neiIdsTrgs)
      }
      else
      {
         for (row in rownames(xcvTrgs))
         {
            d=sqrt(sort(apply(xcvRefs,MARGIN=1,sumSqDiff,xcvTrgs[row,])))[1:object$k]
            neiDstTrgs[row,]=d
            neiIdsTrgs[row,]=names(d)
         }
      }
   }
   else
   {
      for (row in rownames(xTrgs))
      {
         prox=.Call("rfoneprox", INTrefNodes, INTnrow, INTncol,
                     as.integer(nodes[row,]), vector("integer",INTnrow))
         px=sort(prox,index.return = TRUE, decreasing = TRUE)$ix[1:object$k]
         neiDstTrgs[row,]=(INTncol-prox[px])/INTncol
         neiIdsTrgs[row,]=rownames(object$xRefs)[px]
      }
   }
   object$obsDropped=obsDropped
   object$trgRows=trgs
   object$xall=xall
   object$neiDstTrgs=neiDstTrgs
   object$neiIdsTrgs=neiIdsTrgs
   noRefs=TRUE
   object$neiDstRefs=NULL
   object$neiIdsRefs=NULL
   object$ann=ann
   object
}
