// baseline risk

#include <math.h>
#include <stdio.h>
// Weibull
double baselinerisk_weibull(double t, double w1, double w2)
{
  double res;

  if(t == 0)
    res = 0;
  else
    res = pow(w1,2) * pow(w2,2) * pow(t * pow(w1,2), pow(w2,2) - 1); // parametrisation 2
  //res = pow(w1,2) * pow(w2,2) * pow(t, pow(w2,2) - 1); // parametrisation 1
  
  return res;
}


// pour les splines faire lambda(t) = exp(lambda0(t) + Xgamma + Yk)
// et lambda0(t) = sum_l a_l spl_l (t) avec a_l non contraint


// Splines
double baselinerisk_msplines(double t, double * w, double * nodes, int nnodes)
{
  double res;
  int k, l, lm1, lm2, lm3, lp1, lp2, lp3, lp4;
  double ht, htm, h2t, ht2, ht3, hht, h, hh, h2, h3, h4, h3m, h2n, hn, hh2, hh3;

  l = 0;
  for(k = 1; k < nnodes; k++)
    {
      if((t >= nodes[k-1]) & (t < nodes[k]))
	l = k - 1;
    }
  if(t == nodes[nnodes -1])
    l = nnodes - 2;

  res = 0;
  
  if(t == nodes[nnodes - 1])
    {
      res = w[l + 3] * 4.0 / (nodes[l + 1] - nodes[l]);
      //printf("t= %f l=%d wl=%f \n", t, l, w[l]);
    }
  else
    {
      lm1 = l - 1;
      if(l < 2) lm2 = 0; else lm2 = l - 2;  
      if(l < 3) lm3 = 0; else lm3 = l - 3;
      if(l > nnodes - 2) lp1 = nnodes - 1; else lp1 = l + 1;
      if(l > nnodes - 3) lp2 = nnodes - 1; else lp2 = l + 2;
      if(l > nnodes - 4) lp3 = nnodes - 1; else lp3 = l + 3;
      if(l > nnodes - 5) lp4 = nnodes - 1; else lp4 = l + 4;
      
      ht = t - nodes[l];
      htm = t - nodes[lm1];
      h2t = t - nodes[lp2];
      ht2 = nodes[lp1] - t;
      ht3 = nodes[lp3] - t;
      hht = t - nodes[lm2];
      h = nodes[lp1] - nodes[l];
      hh = nodes[lp1] - nodes[lm1];
      h2 = nodes[lp2] - nodes[l];
      h3 = nodes[lp3] - nodes[l];
      h4 = nodes[lp4] - nodes[l];
      h3m = nodes[lp3] - nodes[lm1];
      h2n = nodes[lp2] - nodes[lm1];
      hn = nodes[lp1] - nodes[lm2];
      hh3 = nodes[lp1] - nodes[lm3];
      hh2 = nodes[lp2] - nodes[lm2];

      //printf("t= %f l=%d lp1=%d lp2=%d lp3=%d \n", t, l, lp1, lp2, lp3);
      res = w[l + 3] * (4.0 * ht * ht * ht) / (h4 * h3 * h2 * h) +
	w[l + 2] * ((4.0 * htm * htm * ht2) / (h3m * h2n * hh * h) +
		  (-4.0 * htm * ht * h2t) / (h3m * h2 * h * h2n) +
		  (4.0 * ht3 * ht * ht) / (h3m * h3 * h2 * h)) +
	w[l + 1] * ((4.0 * hht * ht2 * ht2) / (hh2 * hh * h * hn) +
		  (-4.0 * h2t * htm * ht2) / (hh2 * h2n * hh * h) +
		  (4.0 * h2t * h2t * ht) / (hh2 * h2 * h * h2n)) +
	w[l] * (4.0 * ht2 * ht2 * ht2) / (h * hh * hn * hh3);
    }
  
  return exp(res);
  //return res;
}
