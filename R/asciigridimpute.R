AsciiGridPredict = function(object,xfiles,outfiles,xtypes=NULL,lon=NULL,lat=NULL,rows=NULL,cols=NULL,
                            nodata=NULL,myPredFunc=NULL,...)
{
   if (missing(xfiles)   || is.null(xfiles))   stop ("xfiles required")
   if (missing(outfiles) || is.null(outfiles)) stop ("outfiles required")
   if (class(outfiles) != "list") outfiles=list(predict=outfiles)
   if (is.null(names(xfiles))) stop ("xfiles elements must be named")
   if (is.null(names(outfiles))) stop ("outfiles elements must be named")

   return (
      AsciiGridImpute(object,xfiles,outfiles,xtypes=xtypes,ancillaryData=NULL,ann=NULL,
                      lon=lon,lat=lat,rows=rows,cols=cols,nodata=nodata,myPredFunc=myPredFunc,...)
          )
}


AsciiGridImpute = function(object,xfiles,outfiles,xtypes=NULL,ancillaryData=NULL,ann=NULL,
                           lon=NULL,lat=NULL,rows=NULL,cols=NULL,nodata=NULL,myPredFunc=NULL,...)
{
   if (missing(xfiles)   || is.null(xfiles))   stop ("xfiles required")
   if (missing(outfiles) || is.null(outfiles)) stop ("outfiles required")
   if (is.null(names(xfiles))) stop ("xfiles elements must be named")
   if (is.null(names(outfiles))) stop ("outfiles elements must be named")

   factorMatch = get("factorMatch",asNamespace("yaImpute"))

#  make sure there is a type for every xfile

   if (is.null(xtypes))
   {
      xtypes=xfiles
      xtypes[]="numeric"
   }
   else
   {
      tmp=xtypes
      xtypes=xfiles
      xtypes[]=NA
      xtypes[names(tmp)]=tmp
      xtypes[is.na(xtypes)]="numeric"
   }

#  there needs to be an input file for every xvar in the imputation.

   if (is.null(object))
   {
      have=names(xfiles)
   }
   else
   {
      have=intersect(xvars(object),names(xfiles))
      if (length(have) != length(xvars(object)))
      {
         lout = if (length(have)==0) xvars(object) else lout=setdiff(xvars(object),have)
         stop(paste("required maps are missing for variables:",paste(lout ,collapse=", ")))
      }
   }
#  deal with ancillaryData and build allY

   allY = 0
   if (!is.null(ancillaryData) && class(object)=="yai")
   {
      if (length(intersect(class(ancillaryData),c("matrix","data.frame"))) > 0 &&
          nrow(ancillaryData[rownames(object$yRefs),]) == nrow(object$yRefs) &&
          length(intersect(rownames(ancillaryData),rownames(object$yRefs))) == nrow(object$yRefs) )
      {
         toKeep = intersect(union(colnames(ancillaryData),colnames(object$yRefs)),names(outfiles) )
         if (length(toKeep) == 0) allY = NULL
         else
         {
            fromAn=intersect(toKeep,colnames(ancillaryData))
            fromRe=setdiff(intersect(toKeep,colnames(object$yRefs)),fromAn)
            if (length(fromAn)>0 && length(fromRe)>0)
               allY = data.frame(cbind(ancillaryData[rownames(object$yRefs),],object$yRefs)[,toKeep],row.names=rownames(object$yRefs))
            else if (length(fromAn)>0) allY = data.frame(ancillaryData[rownames(object$yRefs),toKeep],row.names=rownames(object$yRefs))
            else if (length(fromRe)>0) allY = data.frame(object$yRefs[,toKeep],row.names=rownames(object$yRefs))
            colnames(allY)=toKeep
         }
      }
      if (is.null(names(allY))) stop ("ancillaryData can not be used")
   }
   if (is.null(names(allY)) && class(object)=="yai")
   {
      toKeep = intersect(colnames(object$yRefs),names(outfiles))
      if (length(toKeep) == 0) allY = NULL
      else allY = data.frame(object$yRefs[,toKeep],row.names=rownames(object$yRefs))
   }
#  if using yai, deal with ann

   if (class(object)=="yai" && is.null(ann)) ann=object$ann

#  set some flags used below

   predYes = length(intersect(names(outfiles),c("predict" ))) == 1
   distYes = length(intersect(names(outfiles),c("distance"))) == 1  && class(object) == "yai"
   useidYes= length(intersect(names(outfiles),c("useid"   ))) == 1  && class(object) == "yai"

   sumIlls=NULL

#  make a list of input file handles and open the files.

   infh = vector("list",length=length(xfiles))
   names(infh)=names(xfiles)
   for (i in 1:length(xfiles))
   {
      infh[[i]]=file(xfiles[[i]])
      open(infh[[i]],open="rt")
   }
   on.exit(lapply(infh,close))

#  make a list of out file handles and open the files.

   if (length(outfiles)> 0)
   {
      outfh = vector("list",length=length(outfiles))
      names(outfh)=names(outfiles)
      for (i in 1:length(outfiles))
      {
         outfh[[i]]=file(outfiles[[i]])
         open(outfh[[i]],open="wt")
      }
      on.exit(lapply(outfh,close),add=TRUE)
   }

#  get and check headers from each input file

   header=NULL
   postWarn=TRUE
   for (i in 1:length(infh))
   {
      newhead = readLines(infh[[i]],n=6)
      if (is.null(header)) header=newhead
      else
      {
         if (!identical(newhead,header))
         {
            cat ("Map headers are not equal\nHeader from file: ",xfiles[[i-1]],"\n")
            print (header)
            cat ("\nHeader from file: ",xfiles[[i]],"\n")
            print (newhead)
            flush.console()
            if (postWarn) warning ("map headers don't match.")
            postWarn=FALSE
         }
         header=newhead
      }
   }

#  write the "common" header to all the output files.

   getVal=function(header,tok)
   {
      for (i in rev(unlist(strsplit(header[grep(tok,header,ignore.case=TRUE)],split=" "))))
      {
         if (i!="") return(as.numeric(i))
      }
   }

   nc  = getVal(header,"NCOLS")
   nr  = getVal(header,"NROWS")
   xllc= getVal(header,"XLLCORNER")
   yllc= getVal(header,"YLLCORNER")
   csz = getVal(header,"CELLSIZE")
   nodv= getVal(header,"NODATA_VALUE")

   if (!is.null(lon)) 
   {
     lon=sort(lon)
     cols=c(max(floor((lon[1]-xllc)/csz),1),min(ceiling((lon[2]-xllc)/csz),nc))
   }
   if (!is.null(lat)) 
   { 
     lat=sort(lat)
     rows=c(max(floor(((yllc+nr*csz)-lat[2])/csz),1),min(nr-ceiling((lat[1]-yllc)/csz),nr))
   }
     
   if (is.null(rows) && is.null(cols) && is.null(nodata)) #header does not change
   {
      for (i in 1:length(outfh))
      {
         cat (header,file=outfh[[i]],sep="\n")
      }
      newnr = nr
      newnc = nc
      nodata= nodv
      rows=c(1,nr)
      cols=c(1,nc)
   }
   else #header changes
   {
      if (is.null(nodata)) nodata=nodv
      if (is.null(rows)) rows=c(1,nr)
      if (is.null(cols)) cols=c(1,nc)
      if (rows[1]>nr) rows[1]=max(1,nr-1)
      if (rows[2]>nr) rows[2]=max(1,nr)
      if (cols[1]>nc) cols[1]=max(1,nc-1)
      if (cols[2]>nc) cols[2]=max(1,nc)
      newnr = rows[2]-rows[1]+1
      newnc = cols[2]-cols[1]+1
      if (rows[1] != 1) yllc = yllc+(csz*rows[1])
      if (cols[1] != 1) xllc = xllc+(csz*cols[1])
      for (i in 1:length(outfh))
      {
         cat("NCOLS         ",as.character(newnc), "\n",file=outfh[[i]],sep="")
         cat("NROWS         ",as.character(newnr), "\n",file=outfh[[i]],sep="")
         cat("XLLCORNER     ",as.character(xllc),  "\n",file=outfh[[i]],sep="")
         cat("YLLCORNER     ",as.character(yllc),  "\n",file=outfh[[i]],sep="")
         cat("CELLSIZE      ",as.character(csz ),  "\n",file=outfh[[i]],sep="")
         cat("NODATA_VALUE  ",as.character(nodata),"\n",file=outfh[[i]],sep="")
      }
   }

   # set up the xlevels. In randomForest version >= 4.5-20, the xlevels
   # are stored in the forest. 
   xlevels = object$xlevels
   if (is.null(xlevels) && class(object) == "randomForest") xlevels=object$forest$xlevels
   if (!is.null(xlevels))
   {  
      if (length(xlevels)>0) for (i in names(xlevels)) if (is.numeric(xlevels[[i]]) && 
                                  length(xlevels[[i]]) == 1) xlevels[[i]] = NULL
      if (length(xlevels) == 0) xlevels=NULL
   }
       
   nskip=0
   if (rows[1]>1) nskip=rows[1]-1

   dpr = max(newnr %/% 100,1)

   cat("Rows per dot: ",dpr," Rows to do:",newnr,"\nToDo: ")

   for (ir in 1:floor(newnr/dpr)) cat (".")
   cat ("\nDone: ");flush.console()

   nodout=suppressWarnings(as.numeric(nodata))

   ircur=0
   for (ir in rows[1]:rows[2])
   {
      indata = vector("list",length=length(xfiles))
      names(indata)=names(infh)
      for (i in 1:length(infh))
      {
         indata[[i]]=scan(infh[[i]],nlines=1,what=vector(mode=xtypes[[i]],length=0),
                     skip=nskip,na.strings=nodv,quiet=TRUE)
      }
      nskip=0
      if (newnc == nc) newdata=data.frame(indata)
      else             newdata=data.frame(indata)[cols[1]:cols[2],,FALSE]
      origRowNames=rownames(newdata)

      if (!is.null(xlevels))
      {
         newdata=factorMatch(newdata,object$xlevels)
         ills = attr(newdata,"illegalLevelCounts")
         if (class(ills)=="list")
         {
            if (is.null(sumIlls))
            {
               sumIlls = ills
               warning ("NA's generated due to illegal level(s).")
            }
            else sumIlls = addIllegalLevels(sumIlls,ills)
         }
      }
      else  attr(newdata,"illegalLevelCounts")=0   # tag the vector so newtargets
                                                   # will not duplicate the data

      newdata=na.omit(newdata)
      omitted=as.vector(attr(newdata,"na.action"))
      if (length(omitted)==length(origRowNames)) # all missing.
      {
         outdata=data.frame(matrix(nodout,length(omitted),length(outfh)),
                            row.names=origRowNames)
         names(outdata)=names(outfh)
      }
      else
      {
         if (!is.null(myPredFunc))
         {
            outdata=myPredFunc(object,newdata,...)
            if (class(outdata) != "data.frame")
               outdata=data.frame(predict=outdata,row.names=rownames(newdata))
         }
         else if (is.null(object))
         {
            outdata=newdata
         }
         else if (class(object) == "yai")
         {
            outdata = NULL
            saveNames=rownames(newdata)
            rownames(newdata)=paste("m",as.character(1:nrow(newdata)),sep="!")
            new = newtargets(object,newdata,ann)
            if (!is.null(allY)) outdata = impute(new,ancillaryData=allY,observed=FALSE)
            rownames(outdata)=saveNames
            if (distYes)  dists = data.frame(distance=new$neiDstTrgs[,1],row.names=rownames(newdata))
            else          dists = NULL
            if (useidYes) useIds= data.frame(distance=new$neiDstTrgs[,1],row.names=rownames(newdata))
            else          useIds= NULL
            if (!is.null(outdata) && !is.null(dists) ) outdata=cbind(outdata,dists)
            else if (is.null(outdata)) outdata=dists
            if (!is.null(outdata) && !is.null(useIds)) outdata=cbind(outdata,useIds)
            else if (is.null(outdata)) outdata=useIds
         }
         else
         {
            predict=predict(object,newdata,...)
            if (is.factor(predict))
            {
               predict = levels(predict)[predict]
               pdnum = as.numeric(predict)
               if (sum(is.na(pdnum)) == 0) predict = pdnum
            }
            outdata=data.frame(predict=predict,row.names=rownames(newdata))
         }
         if (length(omitted)>0)
         {
            # add omitted observations back into data frame in the proper location

            more = data.frame(matrix(nodout,length(omitted),length(names(outdata))),
                              row.names=origRowNames[omitted])
            names(more)=names(outdata)
            outdata = rbind(outdata,more)
            outdata = outdata[sort(as.numeric(rownames(outdata)),index.return = TRUE)$ix,,FALSE]
         }
      }
      for (i in 1:length(outfh))
      {
         vname=names(outfh)[i]
         if (length(intersect(names(outdata),vname))==0) stop (vname," is not present in the prediction")
         if (is.factor(outdata[,vname])) write (as.character(outdata[,vname]),outfh[[i]],ncolumns=newnc)
         else write (outdata[,vname],outfh[[i]],ncolumns=newnc)
      }

      ircur=ircur+1
      if (ircur>=dpr)
      {
         ircur=0
         cat (".");flush.console()
      }
   }
   cat ("\n");flush.console()

   if (!is.null(sumIlls))
   {
      cat ("Summary of illegal levels (those on maps that\nare not present in model fit)\n")
      if (class(sumIlls[[1]]) == "table") print (sumIlls)
      else
      {
         sumIlls = unionDataJoin(sumIlls)
         print (sumIlls)
         cat ("\n")
      }
   }
   sumIlls
}
