#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <map>
using namespace std;

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include "util.h"


void showMatrix(double *x, int xnrow, int xncol){
  int i,j;
  for(i = 0; i < xnrow; i++){
    for(j = 0; j < xncol; j++){
      cout << x[j*xnrow+i] << "\t";
    }
    cout << endl;
  }      
}

SEXP getListElement (SEXP list, const char *str)
{
  SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
  int i;
  
  for (i = 0; i < length(list); i++)
    if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
      elmt = VECTOR_ELT(list, i);
      break;
    }
  
  if(elmt == R_NilValue){
    Rprintf("\nlist element %s not found\n", str);
  }
  return elmt;
}
