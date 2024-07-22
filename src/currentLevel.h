#ifndef CL_H_   
#define CL_H_

double currentLevelGaussian(double * b, double * u, double * x, int nv, int * idef, int * idea);

double currentLevelPoisson(double * b, double * u, double * x, int nv, int * idef, int * idea, int * ido);

#endif
