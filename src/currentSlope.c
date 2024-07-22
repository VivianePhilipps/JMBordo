#include <math.h>
#include "auxfct.h"

// current slope of a marker

// Gaussian
double currentSlopeGaussian(double * b, double * u, double * xdt, int nv, int * idef, int * idea)
{
  double res;

  // slope = mu'(t)
  linpred(b, xdt, u, idef, idea, 1, nv, &res);

  return res;
}

// Poisson
double currentSlopePoisson(double * b, double * u, double * x, double * xdt, int nv, int * idef, int * idea, int * ido)
{
  int k;
  
  // compute mu = Xbeta + Zu
  double mu;
  linpred(b, x, u, idef, idea, 1, nv, &mu);

  // derivative mu'(t)
  double mudt;
  linpred(b, xdt, u, idef, idea, 1, nv, &mudt);

  // slope = mu'(t) * exp(mu(t))
  double res;
  res = mudt * exp(mu);

  // offset
  for(k = 0; k < nv; k++)
    if(ido[k] == 1)
      res = x[k] * res;
  

  return res;
}
