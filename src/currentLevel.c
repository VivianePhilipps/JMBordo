#include <math.h>
#include "auxfct.h"

// current level of a marker

// Gaussian
double currentLevelGaussian(double * b, double * u, double * x, int nv, int * idef, int * idea)
{
  double res;

  // res = Xbeta + Zu
  linpred(b, x, u, idef, idea, 1, nv, &res);

  return res;
}

// Poisson
double currentLevelPoisson(double * b, double * u, double * x, int nv, int * idef, int * idea, int * ido)
{
  int k;
  double res;

  // res = Xbeta + Zu
  linpred(b, x, u, idef, idea, 1, nv, &res);

  // offset
  for(k = 0; k < nv; k++)
    if(ido[k] == 1)
      res = x[k] * res;

  return exp(res);
}
