#ifndef LOGLIK_H_   
#define LOGLIK_H_

// Joint log-likelihoods
void iloglik(double * b, double * bfix, int * fix, int nlong, int nsurv, int * npmlong, int * npmsurv, double * Y, double * Xlong, int * ni, int nvlong, int * idef, int * idea, int * ido, int * typelong, double * Tstart, double * Tstop, int * Event, int * hazard, double * hazardnodes, int * nhazardnodes, double * Xsurv, int * nTi, int nvsurv, int * idsurv, int * asso, int nlinesasso, double * XlongT, double * XlongTdt, int lefftrunc, double * seqQMC, int nQMC, int * indexRE, int * nRE, double * loglik);

//void iloglikJoint(double * b, double * Y, double * X, int nmarkers, int * ni, int nv, int * idef, int * idea, double * seqMC, int nMC, double * loglik);

//void loglikJoint(double * b, double * Y, double * X, int nmarkers, int ns, int * nmes, int nv, int * idef, int * idea, double * seqMC, int nMC, double * loglik);


#endif
