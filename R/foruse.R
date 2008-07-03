# Takes a "yai" object and returns a "for-use" table.
#
# By default (kth=NULL), only the best pick for an observation is listed in the
# use column. For a reference observations, it is always selected to represent
# itself. However, when kth is not NULL, the kth neighbor is reported.
#
# The rowname is a target id, first column is the reference id that
# is select to represent the target, the second is the corresponding
# distance. Every reference is included as a target as well. In those
# cases, the rowname and the use value are the same and the distance
# is zero.

foruse = function (object,kth=NULL)
{
   if (class(object) != "yai") stop ("class must be yai")
   if (!is.null(kth))
   {
      if (kth>object$k) kth=object$k
      if (kth<1)        kth=NULL
   }
   if (!is.null(object$neiIdsRefs))
   {
      if (is.null(kth))
         fu1=data.frame(use=rownames(object$neiIdsRefs),
                        dist=rep(0,nrow(object$neiIdsRefs)),
                        stringsAsFactors = FALSE)
      else
         fu1=data.frame(use=object$neiIdsRefs[,kth],
                        dist=object$neiDstRefs[,kth],
                        stringsAsFactors = FALSE)
      rownames(fu1)=rownames(object$neiIdsRefs)
   }

   else fu1=NULL
   if (!is.null(object$neiIdsTrgs))
   {
      if (is.null(kth)) kth=1
      fu2=data.frame(use=object$neiIdsTrgs[,kth],
                     dist=object$neiDstTrgs[,kth],
                     stringsAsFactors = FALSE)
      rownames(fu2)=rownames(object$neiIdsTrgs)
   }
   else fu2=NULL
   if (is.null(fu1) & is.null(fu2)) return (NULL)
   if      (is.null(fu1)) ans = fu2
   else if (is.null(fu2)) ans = fu1
   else                   ans = rbind(fu1,fu2)
   class(ans)=c("data.frame","foruse.yaImpute")
   ans
}
