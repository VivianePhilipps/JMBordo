#ifndef CS_H_   
#define CS_H_

double currentSlopeGaussian(double * b, double * u, double * x, int nv, int * idef, int * idea);

double currentSlopePoisson(double * b, double * u, double * x, double * xdt,  int nv, int * idef, int * idea, int * ido);

#endif
