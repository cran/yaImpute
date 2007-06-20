//----------------------------------------------------------------------
// File:                rfoneprox.cpp
// Programmer:          Nicholas Crookston
// Description:         computes a proximity vector for one observation
// Last modified:       2007/10/25 (Version 0.2)
//----------------------------------------------------------------------
// Copyright Notice:
// This code was written and prepared by a U.S. Government employee on official 
// time and therefore is in the public domain and not subject to copyright.
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
