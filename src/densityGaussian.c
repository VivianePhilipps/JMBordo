#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "density.h"
#include "auxfct.h"

// density of a Gaussian outcome with random effects u (Monte Carlo points)

// P(Y=y) = 1 / (sigma * sqrt(2 * pi )) * exp( - (y - mu)^2 / (2 * sigma^2))


double densityGaussian(double * b, double * Y, double * X, int ni, int nv, int * idef, int * idea, double * u)
{
  int j, nvc;
  int nef, nea;
  double res, sig2, num;
  double * mu;

  // number of fixed effects
  nef = 0;
  for(j = 0; j < nv; j++)
    {
      if(idef[j] == 1) nef += 1;
    }

  // number of random effects
  nea = 0;
  for(j = 0; j < nv; j++)
    {
      if(idea[j] == 1) nea += 1;
    }
  nvc = 0;
  if(nea > 0) nvc = nea * (nea + 1) / 2;

  // compute mu = Xbeta + Zu
  mu = malloc(ni * sizeof(double));
  linpred(b, X, u, idef, idea, ni, nv, mu);

  // variance of the error
  sig2 = b[nef + nvc] * b[nef + nvc];

  // (Y - mu)^2
  num = 0;
  for (j = 0; j < ni; j++)
    {
      num += (Y[j] - mu[j]) * (Y[j] - mu[j]);
    }

  // compute P(Y=y)
  //loglik = -0.5 * ni * log(2.0*3.14159265358979) - 0.5 * num / sig2 ;
  res = pow(2.0 * 3.14159265358979 * sig2, -ni / 2.0) * exp(-num / (2 * sig2));
   
  free(mu);
  //  printf("densG %f \n",res);
  return res;
}
