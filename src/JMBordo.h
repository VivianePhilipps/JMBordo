#include <R_ext/RS.h>


SEXP loglik_indiv(SEXP b0,
		  SEXP bfix0,
		  SEXP fix0,
		  SEXP nlong0,
		  SEXP nsurv0,
		  SEXP npmlong0,
		  SEXP npmsurv0,
		  SEXP Y0,
		  SEXP Xlong0,
		  SEXP ni0,
		  SEXP nvlong0,
		  SEXP idef0,
		  SEXP idea0,
		  SEXP ido0,
		  SEXP typelong0,
		  SEXP Tstart0,
		  SEXP Tstop0,
		  SEXP Event0,
		  SEXP hazard0,
		  SEXP hazardnodes0,
		  SEXP nhazardnodes0,
		  SEXP Xsurv0,
		  SEXP nTi0,
		  SEXP nvsurv0,
		  SEXP idsurv0,
		  SEXP asso,
		  SEXP nlinesasso,
		  SEXP XlongT0,
		  SEXP XlongTdt0,
		  SEXP lefttrunc0,
		  SEXP seqQMC0,
		  SEXP nQMC0,
		  SEXP indexRE0,
		  SEXP nRE0);
