#include <stdlib.h> 
#include <stdio.h> 

// fonctions auxiliaires

// matrix multiplication
void matmulC (double * m1, double * m2, int l1, int n1, int c2, double * res)
{
  int j, k, l;
  double tmp;

  // m1 : matrice l1*n1 stockee colonne par colonne
  // m2 : matrice n1*c2 stockee colonne par colonne
  // res : matrice l1*c2 stockee colonne par colonne
    
  for (j = 0; j < l1; j++) {
    for (k = 0; k < c2; k++) {
      tmp = 0;
      for (l = 0; l < n1; l++) {
	tmp = tmp + m1[l1 * l + j] * m2[n1 * k + l];
      }
      res[l1 * k + j] = tmp;
    }
  }    

  return;
}

// matrix multiplication
void matmul (double * m1, double * m2, int l1, int n1, int c2, double * res)
{
  int j, k, l;
  double tmp;

  // m1 : matrice l1*n1 stockee ligne par ligne
  // m2 : matrice n1*c2 stockee ligne par ligne
  // res : matrice l1*c2 stockee ligne par ligne
    
  for (j = 0; j < l1; j++) {
    for (k = 0; k < c2; k++) {
      tmp = 0;
      for (l = 0; l < n1; l++) {
	tmp = tmp + m1[n1 * j + l] * m2[c2 * l + k];
      }
      res[c2 * j + k] = tmp;
    }
  }    

  return;
}

// multiplication of a low triangular matrix (mat) and a vector (v)
// result is stored in v
void ltmul(double * mat, double * v, int n)
{
  int j, k;
  double * res = (double *) malloc(n * sizeof(double));

  for(j = 0; j < n; j++)
    {
      res[j] = 0;
      for(k = 0; k <= j; k++)
	res[j] += mat[j * (j + 1) / 2 + k] * v[k];
    }
  
  for(j = 0; j < n; j++)
    v[j] = res[j];

  free(res);

  return;
}


// linear predictor Xbeta + Zu
void linpred(double * b, double * X, double * u, int * idef, int * idea, int ni, int nv, double * res)
{
  int j, k, l, nef, nea;

  // number of fixed effects
  nef = 0;
  for(k = 0; k < nv; k++)
    {
      if(idef[k] == 1)
	{
	  nef += 1;
	}
    }
  
  // number of random effects
  nea = 0;
  for(k = 0; k < nv; k++)
    {
      if(idea[k] == 1)
	{
	  nea += 1;
	}
    }

  // parameters for fixed effects
  double * beta = (double *) malloc(nef * sizeof(double));
  for(k = 0; k < nef; k++)
    {
      beta[k] = b[k];
    }

  // variables in fixed effect
  double * xi = (double *) malloc(nef * ni * sizeof(double));

  for(j = 0; j < ni; j++)
    {
      l = 0;
      for(k = 0; k < nv; k++)
	{
	  if(idef[k] == 1)
	    {
	      xi[j * nef + l] = X[j * nv + k];
	      //printf("k= %d \t j= %d \t X=%f \n",k,j,X[k*nobs + nmescurr*nef + j]);
	      l += 1;
	    }
	}
    }
	  
  // variables in random effect
  double * zi = (double *) malloc(nea * ni * sizeof(double));

  for(j = 0; j < ni; j++)
    {
      l = 0;
      for(k = 0; k < nv; k++)
	{
	  if(idea[k] == 1)
	    {
	      zi[j * nea + l] = X[j * nv + k];
	      l += 1;
	    }
	}
    }

  // res = Xi * beta
  matmul(xi, beta, ni, nef, 1, res);

  // res = Xi * beta + Zi * u
  if(nea > 0)
    {
      double * zu = (double *) malloc(ni * sizeof(double));
      matmul(zi, u, ni, nea, 1, zu);

      for(j = 0; j < ni; j++)
	{
	  res[j] += zu[j];
	}

      free(zu);
    }

  free(beta);
  free(xi);
  free(zi);

  return;
}
