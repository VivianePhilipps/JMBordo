#ifndef BASELINERISK_H_   
#define BASELINERISK_H_

double baselinerisk_weibull(double t, double w1, double w2);

double baselinerisk_msplines(double t, double * w, double * nodes, int nnodes);

#endif
