yai <- function(x=NULL,y=NULL,data=NULL,k=1,noTrgs=FALSE,noRefs=FALSE,
                nVec=NULL,pVal=.05,method="msn",ann=TRUE,mtry=NULL,ntree=500,
                rfMode="buildClasses")
{
   # define functions used internally.
   sumSqDiff=function(x,y) { d=x-y; sum(d*d) }

   findFactors =  get("findFactors",asNamespace("yaImpute"))

   ftest.cor = function (p,q,N,cor)
   {
      s=min(p,q)
      if (s==0) stop ("p and q must be > 0")
      if (length(cor) < s) stop ("cor too short")
      lamda=array(dim=s)
      k=1:s
      for (i in k) lamda[i]=prod(1-cor[i:s]^2)
      r=(N-s-1)-((abs(p-q)+1)/2)
      Ndf=(p-k+1)*(q-k+1)
      u=(Ndf-2)/4
      xx=((p-k+1)^2+(q-k+1)^2)-5
      t=vector(mode="numeric",length=s)
      for (i in k) if (xx[i]>0) t[i]=sqrt(((p-k[i]+1)^2*(q-k[i]+1)^2-4)/xx[i])
      lamda.invt=lamda^(1/t)
      Ddf=(r*t)-(2*u)
      setNA = Ddf < 1 | Ndf < 1
      F=((1-lamda.invt)/lamda.invt)*(Ddf/Ndf) 
      F[setNA] = NA
      pgF=pf(F,Ndf,Ddf,lower.tail=FALSE)
      pgF[setNA] = NA
      if (any(setNA)) warning ("degrees of freedom too low, NA's generated")
      list(F=F,pgF=pgF)
   }
   mymean = function(x)
   {
      if (is.null(ncol(x)))
      {
         ans = if (is.factor(x)) NA else mean(x)
      }
      else
      {
         ans=as.numeric(rep(NA,ncol(x)))
         names(ans)=colnames(x)
         for (i in 1:ncol(x)) if (!is.factor(x[,i])) ans[i]=mean(x[,i])
      }
      ans
   }
   mysd = function(x)
   {
      if (is.null(ncol(x)))
      {
         ans = if (is.factor(x)) NA else sd(x)
      }
      else
      {
         ans=as.numeric(rep(NA,ncol(x)))
         names(ans)=colnames(x)
         for (i in 1:ncol(x)) if (!is.factor(x[,i])) ans[i]=sd(x[,i])
      }
      ans
   }

   #===============================================

   # ARGUMENT and DATA screening

   methodSet=c("msn","msn2","mahalanobis","ica","euclidean","gnn","randomForest","raw",
               "random")
               
   if (!(method %in% methodSet))
      stop (paste("method not one of:",paste(methodSet,collapse =", ")))

   if (method == "gnn") # (GNN), make sure we have package vegan loaded
   {
      if (!require (vegan)) stop("install vegan and try again")
   }
   if (method == "ica") # (ica), make sure we have package fastICA loaded
   {
      if (!require (fastICA)) stop("install fastica and try again")
   }
   if (method == "randomForest") # make sure we have package randomForest loaded
   {
      if (!require (randomForest)) stop("install randomForest and try again")
   }

   cl=match.call()
   obsDropped=NULL
   theFormula=NULL
   yall=NULL
   if (is.data.frame(x) | is.matrix(x))
   {
      if (mode(rownames(x)) != "character") rownames(x)=as.character(rownames(x))
      xall=na.omit (as.data.frame(x))
      if (nrow(xall) != nrow(x))
      {
         warning (nrow(x)-nrow(xall)," x observation(s) removed")
         obsDropped=names(attributes(na.omit(x))$na.action)
      }
      if (!is.null(y))
      {
         if (is.null(dim(y)))
         {
           if (length(y) == nrow (x)) y=data.frame(y,row.names=rownames(x), stringsAsFactors = TRUE)
           else stop("when formulas are not used, y must be a matrix or dataframe, or a vector the same length of rows in x")
         } 
         if (is.matrix(y) | is.data.frame(y))
         {
            if (mode(rownames(y)) != "character") rownames(y)=as.character(rownames(y))
            yall=na.omit(as.data.frame(y))
            if (nrow(yall) != nrow(as.data.frame(y)))
            {
               warning (nrow(y)-nrow(yall)," y observation(s) removed")
               obsDropped=union(obsDropped,names(attributes(na.omit(y))$na.action))
            }
         }
         theFormula=NULL
      }
   }
   else if (class(x) == "formula")
   {
      if (class(y) == "formula") yall=model.frame(y,data=data)
      xall=model.frame(x,data=data)
      obsDropped=setdiff(rownames(data),rownames(xall))
      if (length(obsDropped)) warning (length(obsDropped)," observation(s) removed")
      theFormula=list(x=x,y=y)
   }
   else stop ("x is missing or not a matrix nor dataframe")
   if (is.null(yall) & (method %in% c("mahalanobis","ica","euclidean","randomForest","raw")))
   {
      ydum=TRUE
      yall=data.frame(ydummy=rep(1,nrow(xall)),row.names=rownames(xall))
   }
   else ydum=FALSE
   if (is.null(yall)) stop("y missing")
   if (nrow(xall) == 0) stop ("no observations in x")
   if (! (method %in% c("random","randomForest")))
   {
      fy=0
      if (!(method %in% c("mahalanobis","ica","euclidean","raw"))) fy=sum(findFactors(yall))
      if (fy+sum(findFactors(xall)>0)>0) stop("factors allowed only for methods randomForest or random")
   }
   refs=intersect(rownames(yall),rownames(xall))
   if (length(refs) == 0) stop ("no reference observations.")
   yRefs=yall[refs,,drop=FALSE]
   xRefs=xall[refs,,drop=FALSE]
   trgs=setdiff(rownames(xall),refs)

   if (method == "gnn") # remove rows with zero sums or vegan will error off.
   {
      zero = apply(yRefs,1,sum) <= 0
      ndrop=sum(zero)
      if (ndrop>0)
      {
         warning (ndrop," rows have y-variable row sums <= 0 were converted to target observations for method gnn")
         if (ndrop==length(refs)) stop ("all references were deleted")
         obsDropped=union(obsDropped,refs[zero])
         refs=refs[!zero]
         yRefs=yall[refs,,drop=FALSE]
         xRefs=xall[refs,,drop=FALSE]
         trgs=setdiff(rownames(xall),refs)
      }

      yDrop=apply(yRefs,2,sum) <= 0
      if (sum(yDrop) > 0) warning ("y variables with zero sums: ",
                                    paste(colnames(yRefs)[yDrop],collapse=","))
      if (sum(yDrop) == ncol(yRefs)) stop("no y variables")
      if (sum(yDrop) > 0) yRefs=yRefs[,!yDrop,drop=FALSE]
   }

   # initial scale values (maybe reduced by some methods).
   xScale=list(center=mymean(xRefs),scale=mysd(xRefs))
   yScale=list(center=mymean(yRefs),scale=mysd(yRefs))

   # for all methods except randomForest, random, and raw, variables with zero variance are dropped.
   if (!(method %in% c("randomForest","random","raw")))
   {
      xDrop=xScale$scale < 1e-10
      if (sum(xDrop) > 0) warning ("x variables with zero variance: ",
                                    paste(colnames(xRefs)[xDrop],collapse=","))
      if (sum(xDrop) == ncol(xRefs)) stop("no x variables")
      if (sum(xDrop) > 0)
      {
         xRefs=xRefs[,!xDrop,drop=FALSE]
         xScale$scale=xScale$scale[!xDrop]
         xScale$center=xScale$center[!xDrop]
      }
   }
   else xDrop=NULL

   # for this method, xRefs must be a matrix.
   if (method != "randomForest" && !is.matrix(xRefs)) xRefs=as.matrix(xRefs)

   # define these elements as NULL, some will be redefined below.

   cancor=NULL
   ftest=NULL
   yDrop=NULL
   projector=NULL
   ccaVegan=NULL
   ranForest=NULL
   xTrgs=NULL
   xcvTrgs=NULL
   ICA=NULL

   #======= Define projector (if used), scale the variables, and project the
   # reference space. Also project the target space if it is being used.

   if (method %in% c("msn","msn2")) # msn (both kinds)
   {
      yDrop=yScale$scale < 1e-10
      if (sum(yDrop) > 0) warning ("y variables with zero variance: ",
                                    paste(colnames(yRefs)[yDrop],collapse=","))
      if (sum(yDrop) == ncol(yRefs)) stop("no y variables")
      if (sum(yDrop) > 0)
      {
         yRefs=yRefs[,!yDrop,drop=FALSE]
         yScale$scale=yScale$scale[!yDrop]
         yScale$center=yScale$center[!yDrop]
      }
      xcvRefs=scale(xRefs,center=xScale$center,scale=xScale$scale)
      ycvRefs=scale(yRefs,center=yScale$center,scale=yScale$scale)
      cancor=cancor(xcvRefs,ycvRefs,xcenter = FALSE, ycenter = FALSE)  
                    
      theCols = rownames(cancor$xcoef)

      # scale the coefficients so that the cononical vectors will have unit variance.
      cscal = 1/sd(xcvRefs[,theCols] %*% cancor$xcoef[,1])
      cancor$ycoef = cancor$ycoef * cscal
      cancor$xcoef = cancor$xcoef * cscal
     
      ftest=ftest.cor(p=nrow(cancor$ycoef),q=nrow(cancor$xcoef),N=nrow(yRefs),cancor$cor) 
      if (is.null(nVec)) nVec=length(cancor$cor)-sum(ftest$pgF>pVal)
      if (is.na(nVec)) nVec=1
      nVec=min(nVec,length(cancor$cor))
      nVec=max(nVec,1)
      if (method == "msn" ) projector = cancor$xcoef[,1:nVec,drop=FALSE] %*%
                                        diag(cancor$cor[1:nVec,drop=FALSE],nVec,nVec)
      if (method == "msn2") projector = cancor$xcoef[,1:nVec,drop=FALSE] %*%
                                        diag(cancor$cor[1:nVec,drop=FALSE],nVec,nVec) %*%
                                        diag(sqrt(1/(1-cancor$cor[1:nVec,drop=FALSE]^2)),nVec,nVec)
      if (length(theCols)<ncol(xRefs))
      {
         if (is.null(xDrop)) xDrop=xScale$center==0 #just get the names and create a logical
         remove=setdiff(colnames(xRefs),theCols)
         xDrop[remove]=TRUE
         warning ("x variables with colinearity: ",paste(remove,collapse=","))
         xRefs=xRefs[,theCols,drop=FALSE]
         xScale$center=xScale$center[theCols]
         xScale$scale=xScale$scale[theCols]
         xcvRefs=scale(xRefs,center=xScale$center,scale=xScale$scale)
      }
      xcvRefs=xcvRefs %*% projector
      if (!noTrgs && length(trgs) > 0)
      {
         xTrgs=xall[trgs,theCols,drop=FALSE]
         xcvTrgs=scale(xTrgs,center=xScale$center,scale=xScale$scale)
         xcvTrgs=xcvTrgs %*% projector
      }
   }
   else if (method == "mahalanobis")
   {
      xcvRefs=scale(xRefs,center=xScale$center,scale=xScale$scale)
      qr = qr(xcvRefs)  # maybe we are not at full rank
      xcvRefs=xcvRefs[,qr$pivot[1:qr$rank],drop=FALSE]
      projector = solve(chol(cov(xcvRefs))) # old, wrong, extra transpose removed thanks to Petteri Packalen
      theCols = colnames(projector)
      if (length(theCols)<ncol(xRefs))
      {
         if (is.null(xDrop)) xDrop=xScale$center==0 #just get the names and create a logical
         remove=setdiff(colnames(xRefs),theCols)
         xDrop[remove]=TRUE
         warning ("x variables with colinearity: ",paste(remove,collapse=","))
         xRefs=xRefs[,theCols,drop=FALSE]
         xScale$center=xScale$center[theCols]
         xScale$scale=xScale$scale[theCols]
      }
      nVec = ncol(projector)  # same as qr$rank
      xcvRefs=xcvRefs %*% projector
      if (!noTrgs && length(trgs) > 0)
      {
         xTrgs=xall[trgs,theCols,drop=FALSE]
         xcvTrgs=scale(xTrgs,center=xScale$center,scale=xScale$scale)
         xcvTrgs=xcvTrgs %*% projector
      }
   }
   else if (method == "ica")
   {
      xcvRefs=scale(xRefs,center=xScale$center,scale=xScale$scale)
      qr = qr(xcvRefs)  # maybe we are not at full rank
      xcvRefs=xcvRefs[,qr$pivot[1:qr$rank],drop=FALSE]
      a=fastICA(xcvRefs,ncol(xcvRefs),method="C",)
      ICA=list(S=a$S,K=a$K,A=a$A,W=a$W)
      projector = a$K %*% a$W

      colnames(projector)=colnames(xcvRefs)
      rownames(projector)=colnames(xcvRefs)
      theCols = colnames(xcvRefs)
      if (length(theCols)<ncol(xRefs))
      {
         if (is.null(xDrop)) xDrop=xScale$center==0 #just get the names and create a logical
         remove=setdiff(colnames(xRefs),theCols)
         xDrop[remove]=TRUE
         warning ("x variables with colinearity: ",paste(remove,collapse=","))
         xRefs=xRefs[,theCols,drop=FALSE]
         xScale$center=xScale$center[theCols]
         xScale$scale=xScale$scale[theCols]
      }
      nVec = ncol(projector)  # same as qr$rank
      xcvRefs=xcvRefs %*% projector
      if (!noTrgs && length(trgs) > 0)
      {
         xTrgs=xall[trgs,theCols,drop=FALSE]
         xcvTrgs=scale(xTrgs,center=xScale$center,scale=xScale$scale)
         xcvTrgs=xcvTrgs %*% projector
      }
   }
   else if (method == "euclidean")
   {
      xcvRefs=scale(xRefs,center=xScale$center,scale=xScale$scale)
      nVec = ncol(xRefs)
      if (!noTrgs && length(trgs) > 0)
      {
         xTrgs=xall[trgs,,drop=FALSE]
         xcvTrgs=scale(xTrgs,center=xScale$center,scale=xScale$scale)
      }
   }
   else if (method == "raw")
   {
      xcvRefs=xRefs
      nVec = ncol(xRefs)
      if (!noTrgs && length(trgs) > 0)
      {
         xTrgs=xall[trgs,,drop=FALSE]
         xcvTrgs=as.matrix(xTrgs)
      }
   }
   else if (method == "gnn") # GNN
   {
      xcvRefs=scale(xRefs,center=xScale$center,scale=xScale$scale)
      ccaVegan = cca(X=yRefs, Y=xcvRefs)
      if (is.null(ccaVegan$CCA)) stop ("cca() in package vegan failed, likely cause is too few X or Y variables.")

      # create a projected space for the reference observations
      xcvRefs=predict(ccaVegan,type="lc",rank="full")
      xcvRefs=xcvRefs %*% diag(sqrt(ccaVegan$CCA$eig/sum(ccaVegan$CCA$eig)))

      # create a projected space for the unknowns (target observations)
      if (!noTrgs && length(trgs) > 0)
      {
         xTrgs=xall[trgs,,drop=FALSE]
         xcvTrgs=scale(xTrgs,center=xScale$center,scale=xScale$scale)
         xcvTrgs=predict(ccaVegan,newdata=as.data.frame(xcvTrgs),type="lc",rank="full")
         xcvTrgs=xcvTrgs %*% diag(sqrt(ccaVegan$CCA$eig/sum(ccaVegan$CCA$eig)))
      }
      nVec = ncol(xcvRefs)
   }
   else if (method == "randomForest")
   {  
      rfBuildClasses=NULL
      xTrgs=xall[trgs,1,drop=FALSE]
      rfVersion=packageDescription("randomForest")[["Version"]]  
      if (compareVersion(rfVersion,"4.5-22") < 0) stop("Update your version of randomForest.")
      if (is.null(mtry)) mtry=max(sqrt(ncol(xRefs)),1)
      if (is.null(ntree)) ntree=500
      if (ydum)
      {
         yone=NULL
         ranForest=randomForest(x=xRefs,y=yone,proximity=FALSE,importance=TRUE,keep.forest=TRUE,mtry=mtry,ntree=ntree)
         ranForest$type="yaImputeUnsupervised"
         ranForest=list(unsupervised=ranForest)
      }
      else
      { 
         ranForest=vector("list",ncol(yRefs))
         if (length(ntree) < ncol(yRefs)) ntree=rep(trunc(ntree/ncol(yRefs)),ncol(yRefs))
         for (i in 1:ncol(yRefs))
         {
            yone=yRefs[,i]
            if (!is.factor(yone))
            { 
              if (is.null(rfBuildClasses) && rfMode=="buildClasses") rfBuildClasses=TRUE
              if (is.null(rfBuildClasses))
              {
                 if (compareVersion(rfVersion,"4.5-19") < 0) # if the version is prior to 19
                 {
                   warning("yaImpute directly supports regression for continuous y's for randomForest version 4.5-19 and later.")
                   rfBuildClasses=TRUE
                 }
                 else rfBuildClasses=FALSE
              }
              if (rfBuildClasses)
              {
                yone=as.numeric(yone)
                breaks <- pretty(yone, n = min(20,nclass.Sturges(yone)), min.n = 1)
                div <- diff(breaks)[1]
                yone=as.factor(floor(yone/div))
              }
            }
            ranForest[[i]]=randomForest(x=xRefs,y=yone,proximity=FALSE,importance=TRUE,keep.forest=TRUE,mtry=mtry,ntree=ntree[i])
         }
         names(ranForest)=colnames(yRefs)
      }
      nodes=NULL
      for (i in 1:length(ranForest))
      {
         nodeset=attr(predict(ranForest[[i]],xall,proximity=FALSE,nodes=TRUE),"nodes")
         if (is.null(nodeset)) stop("randomForest did not return nodes")
         colnames(nodeset)=paste(colnames(nodeset),i,sep=".")
         nodes=if (is.null(nodes)) nodeset else cbind(nodes,nodeset)
      }
      refNodes=nodes[rownames(xRefs),]
      INTrefNodes=as.integer(refNodes)
      INTnrow=as.integer(nrow(xRefs))
      INTncol=as.integer(ncol(nodes))
      INTsort = INTrefNodes
      dim(INTsort) = c(INTnrow,INTncol)
      INTsort=apply(INTsort,2,function (x) sort(x,index.return = TRUE, decreasing = FALSE)$ix-1)
      attributes(INTsort)=NULL
      INTsort = as.integer(INTsort)
      attr(ranForest,"rfRefNodeSort") = list(INTrefNodes=INTrefNodes, INTnrow=INTnrow, INTncol=INTncol, INTsort=INTsort)
   }
   else if (method == "random")
   { 
      nVec = 1
      ann=FALSE
      xcvRefs=data.frame(random=runif(nrow(xRefs)),row.names=rownames(xRefs))
      if (!noTrgs && length(trgs) > 0) xcvTrgs=data.frame(random=runif(length(trgs)),row.names=trgs)
   }
   else # default
   {
      stop("no code for specified method")
   }

   k=min(k,nrow(xRefs))

   # ======= find neighbors for TARGETS
   if (noTrgs || length(trgs) == 0)
   {
      neiDstTrgs=NULL
      neiIdsTrgs=NULL
   }
   else
   {
      neiDstTrgs=matrix(data=NA,nrow=length(trgs),ncol=k)
      rownames(neiDstTrgs)=trgs
      colnames(neiDstTrgs)=paste("Dst.k",1:k,sep="")
      neiIdsTrgs=neiDstTrgs
      colnames(neiIdsTrgs)=paste("Id.k",1:k,sep="")
      if (method %in%  c("msn","msn2","mahalanobis","ica","euclidean","gnn","raw"))
      {
         if (ann)
         { 
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
               d=sqrt(sort(apply(xcvRefs,MARGIN=1,sumSqDiff,xcvTrgs[row,]))[1:k])
               neiDstTrgs[row,]=d
               neiIdsTrgs[row,]=names(d)
            }
         }
      }
      else if (method == "random")
      {
         l=k+1
         d = matrix(unlist(lapply(xcvTrgs[[1]],function (x, xcv, l) 
               {
                 sort((xcv-x)^2,index.return=TRUE)$ix[2:l]
               },xcvRefs[[1]],l)),nrow=nrow(xcvTrgs),ncol=k,byrow=TRUE)
         for (ic in 1:ncol(d))
         {
           neiDstTrgs[,ic]=abs(xcvTrgs[,1]-xcvRefs[d[,ic],1])
           neiIdsTrgs[,ic]=rownames(xcvRefs)[d[,ic]]
         }
      }
      else if (method == "randomForest")
      {
        prox=lapply(apply(nodes[rownames(xTrgs),,drop=FALSE],1,as.list),function (x) 
          {
             prx=.Call("rfoneprox", INTrefNodes, INTsort, INTnrow, INTncol,
                       as.integer(x), vector("integer",INTnrow),dup=FALSE) 
             if (k > 1)  px=sort(prx,index.return = TRUE, decreasing = TRUE)$ix[1:k]
             else        px=which.max(prx)
             c(prx[px],px)  # counts followed by pointers to references
          })
        for (i in 1:k)
        {
          neiDstTrgs[,i]=unlist(lapply(prox,function (x,i) (INTncol-x[i])/INTncol,i))
          neiIdsTrgs[,i]=unlist(lapply(prox,function (x,i,k,Rnames) 
                 Rnames[x[k+i]],i,k,rownames(xRefs)))
        } 
      }
      else # default
      {
         stop("no code for specified method")
      }
   }

   # ======= find neighbors for REFERENCES
   if (noRefs)
   {
      neiDstRefs=NULL
      neiIdsRefs=NULL
   }
   else
   {
      neiDstRefs=matrix(data=NA,nrow=nrow(xRefs),ncol=k)
      rownames(neiDstRefs)=rownames(xRefs)
      colnames(neiDstRefs)=paste("Dst.k",1:k,sep="")
      neiIdsRefs=neiDstRefs
      colnames(neiIdsRefs)=paste("Id.k",1:k,sep="")
      l=k+1
      if (method %in%  c("msn","msn2","mahalanobis","ica","euclidean","gnn","raw"))
      {
         if (ann & nrow(xcvRefs)> 0)
         {
             ann.out=ann(xcvRefs, xcvRefs, l, verbose=FALSE)$knnIndexDist
             neiDstRefs[TRUE]=sqrt(ann.out[,(l+2):ncol(ann.out)])
             for (i in 2:l)
                neiIdsRefs[,(i-1)]=rownames(xcvRefs)[ann.out[,i]]
             rownames(neiDstRefs)=rownames(neiIdsRefs)
         }
         else
         {
            for (row in rownames(xcvRefs))
            {
               d=sqrt(sort(apply(xcvRefs,MARGIN=1,sumSqDiff,xcvRefs[row,]))[2:l])
               neiDstRefs[row,]=d
               neiIdsRefs[row,]=names(d)
            }
         }
      }
      else if (method == "randomForest")
      {
        prox=lapply(apply(refNodes,1,as.list),function (x) 
          {
             prx=.Call("rfoneprox", INTrefNodes, INTsort, INTnrow, INTncol,
                       as.integer(x), vector("integer",INTnrow),dup=FALSE) 
             if (k > 1) px=sort(prx,index.return = TRUE, decreasing = TRUE)$ix[2:l]
             else
             { 
               px=which.max(prx)
               prx[px]=-1
               px=which.max(prx)
             }
             c(prx[px],px)  # counts followed by pointers to references
           })
        for (i in 1:k)
        {
           neiDstRefs[,i]=unlist(lapply(prox,function (x,i) (INTncol-x[i])/INTncol,i))
           neiIdsRefs[,i]=unlist(lapply(prox,function (x,i,k,Rnames) 
                  Rnames[x[k+i]],i,k,rownames(xRefs)))
        } 
      }
      else if (method == "random")
      {
         l=k+1
         d = matrix(unlist(lapply(xcvRefs[[1]],function (x, xcv, l) 
               {
                 sort((xcv-x)^2,index.return=TRUE)$ix[2:l]
               },xcvRefs[[1]],l)),nrow=nrow(xcvRefs),ncol=k,byrow=TRUE)
               
         for (ic in 1:ncol(d))
         {
           neiDstRefs[,ic]=abs(xcvRefs[,1]-xcvRefs[d[,ic],1])
           neiIdsRefs[,ic]=rownames(xcvRefs)[d[,ic]]
         }
      }
      else # default
      {
         stop("no code for specified method")
      }
   }

   xlevels=NULL
   fa=findFactors(xRefs)
   if (sum(fa)>0)
   {
      xlevels=vector(mode="list",length=sum(fa))
      k=0
      for (i in 1:length(fa))
      {
         if (fa[i])
         {
            k=k+1
            xlevels[[k]]=levels(xRefs[,i])
            names(xlevels)[[k]]=names(xRefs)[i]
         }
      }
   }

   out=list(call=cl,yRefs=yRefs,xRefs=xRefs,obsDropped=obsDropped,yDrop=yDrop,
            xDrop=xDrop,trgRows=trgs,xall=xall,cancor=cancor,theFormula=theFormula,
            ftest=ftest,yScale=yScale,xScale=xScale,ccaVegan=ccaVegan,ranForest=ranForest,
            ICA=ICA,k=k,projector=projector,nVec=nVec,pVal=pVal,method=method,ann=ann,
            xlevels=xlevels,neiDstTrgs=neiDstTrgs,neiIdsTrgs=neiIdsTrgs,
            neiDstRefs=neiDstRefs,neiIdsRefs=neiIdsRefs)

   class(out)="yai"
   out
}
