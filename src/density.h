#ifndef DENSITY_H_   
#define DENSITY_H_


// Univariate distributions
double densityGaussian(double * b, double * Y, double * X, int ni, int nv, int * idef, int * idea, double * u);

double densityPoisson(double * b, int * Y, double * X, int ni, int nv, int * idef, int * idea, int * ido, double * u);


#endif
