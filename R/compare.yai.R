compare.yai = function (...,ancillaryData=NULL,vars=NULL,method="rmsd")
{
   if (missing(...)) stop ("... required")
   okClasses <- c("yai","impute.yai")
   for (object in list(...))
      if (length(intersect(class(object),okClasses))==0)
         stop("object classes must be one of ",paste(okClasses,collapse=", "))
   args <- list(...)
   if (length(intersect(method,c("rmsd","cor")))==0) stop("method must be rmsd or cor")
   names(args) <- as.list(substitute(list(...)))[-1]  #who would of guessed that this is the way!
   ans <- vector("list",length(args))
   i <- 0
   tag <- if (method=="rmsd") "rmsdS" else "cor"

   for (object in list(...))
   {
      i <- i+1
      if (inherits(object,"yai")) object <- impute.yai(object,ancillaryData=ancillaryData,vars=vars,observed=TRUE)
      one <- switch(match(method,c("rmsd","cor")),
                    rmsd.yai(object,vars=vars,scale=TRUE),
                    cor.yai(object,vars=vars), NULL)
      names(one) <- paste(names(args)[i],tag,sep=".")
      ans[[i]] <- one
   }
   names(ans) <- names(args)
   ans <- unionDataJoin(ans)
   class(ans) <- c("compare.yai",class(ans))
   ans
}
