#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "JMBordo.h"

static R_CallMethodDef CallEntries[] = {
  {"loglik_indiv", (DL_FUNC) &loglik_indiv, 34},
  {NULL, NULL, 0}
};



void R_init_JMBordo(DllInfo * dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
