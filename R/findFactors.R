findFactors = function(x)
{
   if (is.null(x)) return(NULL)
   if (is.data.frame(x))
   {
      factors=rep(FALSE,ncol(x))
      for (i in 1:ncol(x)) factors[i]=is.factor(x[,i])
   }
   else if (is.matrix(x)) factors=rep(FALSE,ncol(x))
   else factors=is.factor(x)
   factors
}

factorMatch = function (x,xlevels)
{
   if (is.matrix(x) || !is.data.frame(x)) return(x = x)

   if (!is.null(attr(x,"illegalLevelCounts"))) return(x = x)
   facts=intersect(names(x),names(xlevels))
   if (length(facts)==0) return(x)
   nas=NULL
   miss=NULL
   for (varName in facts)
   {
      origLevs= unique(I(as.character(x[,varName])))
      new = match(x[,varName],xlevels[[varName]])
      nas=is.na(new)
      if (sum(nas)>0)
      {
         class(new)="factor"
         attr(new,"levels")=xlevels[[varName]]
         mtb=table(as.character(x[nas,varName]))
         miss=if(is.null(miss)) list(mtb) else c(miss,list(mtb))
         x[,varName]=new
      }
   }
   names(miss)=facts
   attr(x,"illegalLevelCounts")=if (is.null(miss)) 0 else miss
   x
}

addIllegalLevels=function(a1,a2)
{
   if (class(a1)=="data.frame") a1=attr(a1,"illegalLevelCounts")
   if (class(a2)=="data.frame") a2=attr(a2,"illegalLevelCounts")
   vars=union(names(a1),names(a2))
   if (length(vars)==0) return(NULL)
   out =vector(mode = "list", length = length(vars))
   names(out)=vars
   for (var in vars)
   {
      m1 = a1[[var]]
      if (!is.null(m1))
      {
         m1 = as.matrix(m1)
         colnames(m1)="m1"
      }
      m2 = a2[[var]]
      if (!is.null(m2))
      {
         m2 = as.matrix(m2)
         colnames(m2)="m2"
      }
      if (is.null(m1) && is.null(m2)) return (NULL)
      else if (is.null(m1)) both = m2
      else if (is.null(m2)) both = m1
      else both = unionDataJoin(m1,m2)
      both[is.na(both)]=0
      both = as.matrix(apply (both,1,sum))
      colnames(both)=var
      out[[var]]=both
   }
   out
}

