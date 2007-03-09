unionDataJoin=function(...,warn=TRUE)
{
#  creates a data frame that has the rows defined by a union of all rownames in the
#  arguments and columns defined by a union of all colnames in the arguments.
#  a single argument can be a list of data frames or matrices
#
#  when warn is TRUE, columns that occur in more than one source are listed in a warning.

   args=list(...)
   if (length(args)==1)
   {
      args=args[[1]]
      if (class(args) != "list" ) stop("single argument must be a list")
   }
   for (d in args)
   {
      if (!is.data.frame(d) & !is.matrix(d)) stop ("arguments or list members must be matrices or data frames")
      if (is.matrix(d))
     {
         if (is.null(colnames(d))) stop ("column names are requried within all input matrices")
         if (is.null(rownames(d))) stop ("row names are requried within all input matrices")
         if (length(unique(colnames(d))) != length(colnames(d))) stop("column names must be unique within all input matrices")
      }
   }
   rows=NULL
   cols=NULL
   haveCol=NULL
   for (d in args)
   {
      rows=union(rows,rownames(d))
      haveCol=union(intersect(cols,colnames(d)),haveCol)
      cols=union(cols,colnames(d))
   }
   if (warn & length(haveCol)>0)
      warning ("Columns: \"",paste(haveCol,collapse=", "),
               "\" were defined more than once")
   all=as.data.frame(matrix(data=NA,nrow=length(rows),ncol=length(cols)))
   rownames(all)=rows
   colnames(all)=cols
   factors=matrix(data=FALSE,nrow=length(cols),ncol=1)
   rownames(factors)=cols
   colnames(factors)="factor"
   for (d in args)
   {
      theCols=colnames(d)
      both = intersect(rownames(all),rownames(d))
      if (is.data.frame(d))
      {
         for (var in theCols)
         {
            if (is.factor(d[,var]))
            {
               factors[var,1] = TRUE
               all[both,var]=levels(d[,var])[d[both,var]]
            }
            else all[both,var]=d[both,var]
         }
      }
      else all[both,theCols]=d[both,]
   }
   for (var in colnames(all)) if (factors[var,1]) all[,var]=as.factor(all[,var])
   all
}