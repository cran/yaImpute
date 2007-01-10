//----------------------------------------------------------------------
// File:                        rfoneprox.c
// Programmer:          Nicholas Crookston
// Description:         computes a proximity vector for one observation
// Last modified:       2006/06/30 (Version 0.1)
//----------------------------------------------------------------------
// Copyright (c) None
//
//----------------------------------------------------------------------
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

extern "C" {
 SEXP rfoneprox(SEXP nodes, SEXP nobs, SEXP ntree, SEXP obs, SEXP prox)
 {
   int itr, myNode, iob;
   for (itr=0; itr<INTEGER(ntree)[0]; itr++)
   {
      myNode=INTEGER(obs)[itr];
      for (iob=0; iob<INTEGER(nobs)[0]; iob++)
      {
         if (INTEGER(nodes)[INTEGER(nobs)[0]*itr+iob] == myNode) 
		 INTEGER(prox)[iob]++;
      }
   }
 return(prox);
 }
}
