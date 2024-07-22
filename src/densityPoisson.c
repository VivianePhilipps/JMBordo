#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "density.h"
#include "auxfct.h"

// density of a Poisson outcome with random effects u (Monte Carlo points)

// P(Y = y) = (eta^y * exp(-eta)) / y!
// and log(eta) = Xbeta + Zu 
// or log(eta) = log(offset) + Xbeta + Zu


double densityPoisson(double * b, int * Y, double * X, int ni, int nv, int * idef, int * idea, int * ido, double * u)
{
  int j, k;
  int offset;
  double res, fact;
  double * offsetvalue;
  
  // compute mu = Xbeta + Zu
  double * mu = (double *) malloc(ni * sizeof(double));
  linpred(b, X, u, idef, idea, ni, nv, mu);

  // check is offset
  offset = 0;
  for(k = 0; k < nv; k++)
    {
      if(ido[k] == 1)
	{
	  offset = 1;

	  offsetvalue = (double *) malloc(ni * sizeof(double));
	  for(j = 0; j < ni; j++)
	    offsetvalue[j] = X[nv * k + j];

	  break;
	}
    }

  // compute log(P(Y=y))
  res = 1;
  for(j = 0; j < ni; j++)
    {
      // y!
      fact = 1;
      for(k = 0; k < Y[j]; k++)
	fact *= k;
      res = res / fact;

      // eta^y = exp(mu)^y = exp(mu*y)
      res *= exp(mu[j] * Y[j]);

      // exp(-eta) = exp(-exp(mu))
      if(offset == 0)
	res *= exp(-exp(mu[j]));
      else
	res *= pow(offsetvalue[j], Y[j]) * exp(-offsetvalue[j] * exp(mu[j]));
    }

  if(offset == 1) free(offsetvalue);
  free(mu);
  
  return res;
}
