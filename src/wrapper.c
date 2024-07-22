#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include "loglik.h"


SEXP loglik_indiv(SEXP b0, SEXP bfix0, SEXP fix0, SEXP nlong0, SEXP nsurv0, SEXP npmlong0, SEXP npmsurv0,  SEXP Y0, SEXP Xlong0, SEXP ni0, SEXP nvlong0, SEXP idef0, SEXP idea0, SEXP ido0, SEXP typelong0, SEXP Tstart0, SEXP Tstop0, SEXP Event0, SEXP hazard0, SEXP hazardnodes0, SEXP nhazardnodes0, SEXP Xsurv0, SEXP nTi0, SEXP nvsurv0, SEXP idsurv0, SEXP asso0, SEXP nlinesasso0, SEXP XlongT0, SEXP XlongTdt0, SEXP lefttrunc0, SEXP seqQMC0, SEXP nQMC0, SEXP indexRE0, SEXP nRE0)
{
  double * b = REAL(b0); // REAL et INTEGER renvoient des pointeurs
  double * bfix = REAL(bfix0);
  int * fix = INTEGER(fix0);
  int nlong = INTEGER(nlong0)[0];
  int nsurv = INTEGER(nsurv0)[0];
  int * npmlong = INTEGER(npmlong0);
  int * npmsurv = INTEGER(npmsurv0);
  double * Y = REAL(Y0);
  double * Xlong = REAL(Xlong0);
  int * ni = INTEGER(ni0);
  int nvlong = INTEGER(nvlong0)[0];
  int * idef = INTEGER(idef0);
  int * idea = INTEGER(idea0);
  int * ido = INTEGER(ido0);
  int * typelong = INTEGER(typelong0);
  double * Tstart = REAL(Tstart0);
  double * Tstop = REAL(Tstop0);
  int * Event = INTEGER(Event0);
  int * hazard = INTEGER(hazard0);
  double * hazardnodes = REAL(hazardnodes0);
  int * nhazardnodes = INTEGER(nhazardnodes0);
  double * Xsurv = REAL(Xsurv0);
  int * nTi = INTEGER(nTi0);
  int nvsurv = INTEGER(nvsurv0)[0];
  int * idsurv = INTEGER(idsurv0);
  int * asso = INTEGER(asso0);
  int nlinesasso = INTEGER(nlinesasso0)[0];
  double * XlongT = REAL(XlongT0);
  double * XlongTdt = REAL(XlongTdt0);
  int lefttrunc = INTEGER(lefttrunc0)[0];
  double * seqQMC = REAL(seqQMC0);
  int nQMC = INTEGER(nQMC0)[0];
  int * indexRE = INTEGER(indexRE0);
  int * nRE = INTEGER(nRE0);

  double loglik = 0;
  iloglik(b, bfix, fix, nlong, nsurv, npmlong, npmsurv, Y, Xlong, ni, nvlong, idef, idea, ido, typelong, Tstart, Tstop, Event, hazard, hazardnodes, nhazardnodes, Xsurv, nTi, nvsurv, idsurv, asso, nlinesasso, XlongT, XlongTdt, lefttrunc, seqQMC, nQMC, indexRE, nRE, &loglik);
  //  printf("dans wrap, loglik=%f \n",loglik);
  SEXP loglik0;
  PROTECT(loglik0=allocVector(REALSXP,1));
  REAL(loglik0)[0] = loglik;

  UNPROTECT(1);
  return loglik0;
}



/* SEXP loglik_Joint (SEXP b0, SEXP Y0, SEXP X0, SEXP nmarkers0, SEXP ns0, SEXP nmes0, SEXP nv0, SEXP idef0, SEXP idea0, SEXP seqMC0, SEXP nMC0) */
/* { */
/*   double * b = REAL(b0); */
/*   double * Y = REAL(Y0); */
/*   double * X = REAL(X0); */
/*   int nmarkers = INTEGER(nmarkers0)[0]; */
/*   int ns = INTEGER(ns0)[0]; */
/*   int * nmes = INTEGER(nmes0)[0]; */
/*   int nv = INTEGER(nv0)[0]; */
/*   int * idef = INTEGER(idef0)[0]; */
/*   int * idea = INTEGER(idea0)[0]; */
/*   int * seqMC = INTEGER(seqMC0)[0]; */
/*   int nMC = INTEGER(nMC0)[0]; */

/*   double loglik = 0; */
/*   loglikJoint(b, Y, X, nmarkers, nmarkers, ns, nmes, nv, idef, idea, seqMC, nMC, &loglik); */
/*   //  printf("dans wrap, loglik=%f \n",loglik); */
/*   SEXP loglik0; */
/*   PROTECT(loglik0=allocVector(REALSXP,1)); */
/*   REAL(loglik0)[0] = loglik; */

/*   UNPROTECT(1); */
/*   return loglik0; */
/* } */




