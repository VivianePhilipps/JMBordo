#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "loglik.h"
#include "auxfct.h"
#include "density.h"
#include "baselinerisk.h"
#include "currentLevel.h"
#include "currentSlope.h"
#include "nonlinear.h"

// Xlong ligne par ligne
// plusieurs lignes surv par sujet. Voir si nTi=0 ok.
// on supp que asso est trie par evt et par Y


void iloglik(double * b, double * bfix, int * fix, int nlong, int nsurv, int * npmlong, int * npmsurv, double * Y, double * Xlong, int * ni, int nvlong, int * idef, int * idea, int * ido, int * typelong, double * Tstart, double * Tstop, int * Event, int * hazard, double * hazardnodes, int * nhazardnodes, double * Xsurv, int * nTi, int nvsurv, int * idsurv, int * asso, int nlinesasso, double * XlongT, double * XlongTdt, int lefttrunc, double * seqQMC, int nQMC, int * indexRE, int * nRE, double * loglik)
{
  int l, sumnpm1, sumnpm2, k, j, jfix, sumQMC, sumni, sumnTi, i, m, e, ke, p;
  int nes, nrisq, sumnRE, np, avtlong, ncovRE, sumasso, sumhazardnodes, jasso, jasso2;
  double vraisl, vraisk;
  int npmtot, nREtot;
  double T0, T, fevt, xegamma, baserisk, cumhaz, t;
  int D;
  double yt, ytdt;
  double alphakasso, alphaAsso;
  double trunccumhaz, truncl, trunc;
  int * idefk;
  int * ideak;
  int * idok;
  int * idsurve;
  double * nodes;
  double * xe;
  double * be;
  double * bk;
  double * uk;
  double * xkT;
  double * xkTdt;
  double * yk;
  int * intyk;
  double * xk;
  double * ptsGK;
  double * wGK;

  // total number of random effects (RE)
  nREtot = 0;
  for(k = 0; k < nlong; k++)
    nREtot += nRE[k];

  // covariance between markers
  ncovRE = 0;
  if(nREtot > 0) ncovRE = nREtot * (nREtot + 1 ) /2;
  
  // total number of parameters (estimated + fixed)
  npmtot = 0;
  avtlong = 0;
  if(nlong > 0)
    {
      for(k = 0; k < nlong; k++)
	{
	  npmtot += npmlong[k];
	  if(nRE[k] > 0) ncovRE -= nRE[k] * (nRE[k] + 1) / 2;
	}
    }
  if(nsurv > 0)
    {
      for(k = 0; k < nsurv; k++)
	{
	  npmtot += npmsurv[k];
	  avtlong += npmsurv[k];
	}
    }
  npmtot += ncovRE;

  // total vector of parameters (b and bfix)
  double * btot = (double *) malloc(npmtot * sizeof(double));
  // double * btot =  malloc(npmtot * sizeof * btot); mieux? !**
  j = 0;
  jfix = 0;
  for(k = 0; k < npmtot; k++)
    {
      if(fix[k] == 0)
	{
	  btot[k] = b[j];
	  j += 1;
	}
      else
	{
	  btot[k] = bfix[jfix]; 
	  jfix += 1;
	}
    }
  
  // cholesky (lower triangular) of the variance matrix of RE
  double * bRE = malloc(nREtot * (nREtot + 1) / 2 * sizeof(double)); 
  for(k = 0; k < nREtot * (nREtot + 1) / 2; k++)
    {
      bRE[k] = btot[indexRE[k]];
    }

  // allocation of U
  double * U = (double *) malloc(nREtot * sizeof(double));

  // Gauss-Kronrod
  if(nsurv > 0)
    {
      ptsGK = (double *) malloc(15 * sizeof(double));
      wGK = (double *) malloc(15 * sizeof(double));

      ptsGK[0] = 0;
      ptsGK[1] = 0.207784955007898467600689403773245;
      ptsGK[2] = 0.405845151377397166906606412076961;
      ptsGK[3] = 0.586087235467691130294144838258730;
      ptsGK[4] = 0.741531185599394439863864773280788;
      ptsGK[5] = 0.864864423359769072789712788640926;
      ptsGK[6] = 0.949107912342758524526189684047851;
      ptsGK[7] = 0.991455371120812639206854697526329;
      ptsGK[8] = -0.207784955007898467600689403773245;      
      ptsGK[9] = -0.405845151377397166906606412076961;
      ptsGK[10] = -0.586087235467691130294144838258730;
      ptsGK[11] = -0.741531185599394439863864773280788;
      ptsGK[12] = -0.864864423359769072789712788640926;
      ptsGK[13] = -0.949107912342758524526189684047851;
      ptsGK[14] = -0.991455371120812639206854697526329;

      wGK[0] = 0.209482141084727828012999174891714;
      wGK[1] = 0.204432940075298892414161999234649;
      wGK[2] = 0.190350578064785409913256402421014;
      wGK[3] = 0.169004726639267902826583426598550;
      wGK[4] = 0.140653259715525918745189590510238;
      wGK[5] = 0.104790010322250183839876322541518;
      wGK[6] = 0.063092092629978553290700663189204;
      wGK[7] = 0.022935322010529224963732008058970;
      wGK[8] = 0.204432940075298892414161999234649;
      wGK[9] = 0.190350578064785409913256402421014;
      wGK[10] = 0.169004726639267902826583426598550;
      wGK[11] = 0.140653259715525918745189590510238;
      wGK[12] = 0.104790010322250183839876322541518;
      wGK[13] = 0.063092092629978553290700663189204;
      wGK[14] = 0.022935322010529224963732008058970;	  
    }


  // initialize result to 0
  *loglik = 0;
  trunc = 0;
  
  // QMC loop for integration over the RE ( loglik = sum_l f(Y | RE_l) * f(T,D | RE_l) )
  sumQMC = 0;
  for(l = 0; l < nQMC; l++)
    {
      vraisl = 1.0;
      truncl = 1.0;
      
      // RE for iteration l
      if(nREtot > 0)
	{
	  for(k = 0; k < nREtot; k++)
	    U[k] = seqQMC[sumQMC + k]; // in entry U is from a standard Gaussian
	  
	  ltmul(bRE, U, nREtot); // here U has the right variance 
	}

      //   if(l < 2){
      //printf("bRE = \n ");
      //for(k=0; k<nREtot; k++) printf("%f et U %f\t",bRE[k],U[k]);
      //printf("\n ");
      //}
      // survival f(T,D | RE_l)
      if(nsurv > 0)
	{
	  // avec un produit sur les evts/lignes
	  // donc Xsurv avec plusieurs lignes par sujets/
	  

	  idsurve = (int *) malloc(nvsurv * sizeof(int));
	  xe = (double *) malloc(nvsurv * sizeof(double));

	  fevt = 1.0;

	  sumnpm1 = 0;
	  sumhazardnodes = 0;
	  sumnTi = 0;
	  jasso = 0;
	  jasso2 = 0;
	  for(e = 0; e < nsurv; e++)
	    {
	      // parameters for event e
	      be = (double *) malloc(npmsurv[e] * sizeof(double));
	      for(i = 0; i < npmsurv[e]; i++)
		be[i] = btot[sumnpm1 + i];
	      
	      // nodes if splines
	      if(hazard[e] == 2)
		{
		  nodes = (double *) malloc(nhazardnodes[e] * sizeof(double));
		  
		  for(m = 0; m < nhazardnodes[e]; m++)
		    nodes[m] = hazardnodes[sumhazardnodes + m];
		}
	      
	      for(ke = 0; ke < nTi[e]; ke++)
		{
		  // outcomes for event e
		  T0 = Tstart[sumnTi + ke];
		  T = Tstop[sumnTi + ke];
		  D = Event[sumnTi + ke];

		  /* // skip the computations if D is -1 */
		  /* if (D == -1) */
		  /* 	{ */
		  /* 	  sumnpm1 += npmsurv[e]; */
		  /* 	  continue; //? pour le multi-etat, transition pas possible */
		  /* 	} */	      
	      
		  // covariates for event e
		  for(m = 0; m < nvsurv; m++)
		    {
		      xe[m] = Xsurv[nvsurv * sumnTi + m];
		      idsurve[m] = idsurv[nvsurv * sumnTi + m];
		    }
		  
		  
		  // hazard at time t
		  cumhaz = 0;
		  trunccumhaz = 0;
		  np = 16;
		  if(lefttrunc == 1) np += 15;
		  for(p = 0; p < np; p++)
		    {
		      // skip computation if no event
		      if((p == 0) & (D == 0)) continue;
		      
		      // time t
		      if(p == 0) 
			t = T; //  event time
		      else if((p > 0) & (p < 16))
			t = (ptsGK[p - 1] + 1) / 2 * T; // integration point 0->T
		      else
			t = (ptsGK[p - 16] + 1) / 2 * T0;// integration point 0->T0
		      
		      // baseline risk
		      nrisq = 0;
		      baserisk = 0;
		      if(hazard[e] == 1)
			{
			  nrisq = 2;
			  baserisk = baselinerisk_weibull(t, be[0], be[1]);
			}
		      else if(hazard[e] == 2)
			{
			  nrisq = nhazardnodes[e] + 2;
			  baserisk = baselinerisk_msplines(t, be, nodes, nhazardnodes[e]);
			}
		      
      /* if(l < 2){ */
      /* 	printf("baselinerisk = %f en t=%f \n ", baserisk, t); */
      /* printf("\n "); */
      /* } */
		      // time fixed covariates
		      if(p <= 1)
			{
			  nes = 0; // number of covariates in survival, ie sum(idsurve)
			  xegamma = 0; // xe * gamma
			  for(m = 0; m < nvsurv; m++)
			    {
			      if(idsurve[m] == 1)
				{
				  xegamma += xe[m] * be[nrisq + nes];
				  nes++;
				}
			    }
			  /* if(l < 2){ */
			  /* 	printf("xegamma = %f \n ", xegamma); */
			  /* printf("\n "); */
			  /* } */
			}
		  
		      // longitudinal covariate in survival
		      alphaAsso = 0;
		      yt = 0 ;
		      ytdt = 0;
		      jasso = jasso2;
		      if(nlong > 0)
			{
			  // allocation of indicators of fixed and random effects for marker k
			  idefk = (int *) malloc(nvlong * sizeof(int));
			  ideak = (int *) malloc(nvlong * sizeof(int));
			  idok = (int *) malloc(nvlong * sizeof(int));
			  
			  // allocation of vector of covariates
			  xkT = (double *) malloc(nvlong * sizeof(double));
			  
			  sumnRE = 0;
			  sumnpm2 = avtlong;
			  sumasso = 0;
			  alphakasso = 0;
			  for(k = 0; k < nlong; k++)
			    {
			      // check if marker k is associated to event e
			      if((asso[5 * jasso + 1] != k+1) | (asso[5 * jasso + 2] != e+1))
				{
				  sumnpm2 += npmlong[k];
				  sumnRE += nRE[k];
				  continue;
				}
			    
			      // parameters for marker k
			      bk = (double *) malloc(npmlong[k] * sizeof(double));
			      for(i = 0; i < npmlong[k]; i++)
				bk[i] = btot[sumnpm2 + i];
			      
			      // indicator of fixed and random effects for marker k
			      for(j = 0; j < nvlong; j++)
				{
				  idefk[j] = idef[k * nvlong + j];
				  ideak[j] = idea[k * nvlong + j];
				  idok[j] = ido[k * nvlong + j];
				}
			  
			      // RE for marker k
			      uk = (double *) malloc(nRE[k] * sizeof(double));
			      for(i = 0; i < nRE[k]; i++)
				uk[i] = U[sumnRE + i];
			      
			  
			      // Yk's associations with risk function at time t :
			      alphakasso = 0;
			      while((jasso < nlinesasso) & (asso[5 * jasso + 1] == k+1) & (asso[5 * jasso + 2] == e+1))
				{
				  // random effects
				  if(asso[5 * jasso] == 1)
				    {
				      if(asso[5 * jasso + 3] == 1) //linear effect
					{
					  for(i = 0; i < nRE[k]; i++)
					    alphakasso += uk[i] * be[nrisq + nes + sumasso + i];
					  
					  sumasso += nRE[k];
					}
				      else if(asso[5 * jasso + 3] == 2) //polynomial nonlinear effect
					{
					  for(i = 0; i < nRE[k]; i++)
					    {
					      alphakasso += polynomials(uk[i], be[nrisq + nes + sumasso], be[nrisq + nes + sumasso + 1], be[nrisq + nes + sumasso + 2]);
					      
					      sumasso += 3;
					    }
					}
				      else if(asso[5 * jasso + 3] == 3) //logistic nonlinear effect
					{
					  for(i = 0; i < nRE[k]; i++)
					    {				      
					      alphakasso = logistic(uk[i], be[nrisq + nes + sumasso], be[nrisq + nes + sumasso + 1], be[nrisq + nes + sumasso + 2]);
					      
					      sumasso += 3;
					    }				      
					}
				    }
				 
				  // current Level
				  if(asso[5 * jasso] == 2)
				    {
				      // covariates from XlongT
				      for(m = 0; m < nvlong; m++)
					{
					  //xkT[m] = XlongT[nvlong * (np * (nlong * e + k) + p)  + m];
					  xkT[m] = XlongT[nvlong * (sumnTi * np * nlong + k * np + p)  + m];
					}
				      
				      // current value of marker k
				      if(typelong[k] == 1)
					{
					  yt = currentLevelGaussian(bk, uk, xkT, nvlong, idefk, ideak);
					}
				      else if(typelong[k] == 2)
					{
					  yt = currentLevelPoisson(bk, uk, xkT, nvlong, idefk, ideak, idok);
					}
				      else if(typelong[k] == 3)
					{
					  //add...
					}
				 
				      if(asso[5 * jasso + 3] == 1) //linear effect
					{
					  
					  alphakasso = be[nrisq + nes + sumasso] * yt;
					  
					  sumasso += 1;
					}
				      else if(asso[5 * jasso + 3] == 2) //polynomial nonlinear effect
					{
					  alphakasso = polynomials(yt, be[nrisq + nes + sumasso], be[nrisq + nes + sumasso + 1], be[nrisq + nes + sumasso + 2]);
					  sumasso += 3;
					}
				      else if(asso[5 * jasso + 3] == 3) //logitic nonlinear effect
					{
					  alphakasso = logistic(yt, be[nrisq + nes + sumasso], be[nrisq + nes + sumasso + 1], be[nrisq + nes + sumasso + 2]);
					  
					  sumasso += 3;
					}
				    }
				  
				  // current Slope
				  if(asso[5 * jasso] == 3)
				    {
				      // covariates from XlongTdt
				      xkTdt = (double *) malloc(nvlong * sizeof(double));
				      for(m = 0; m < nvlong; m++)
					{
					  xkTdt[m] = XlongTdt[nvlong * (sumnTi * np * nlong + k * np + p)  + m];
					}
			      
				      // current slope of marker k
				      if(typelong[k] == 1)
					{
					  ytdt = currentSlopeGaussian(bk, uk, xkTdt, nvlong, idefk, ideak);
					}
				      else if(typelong[k] == 2)
					{
					  ytdt = currentSlopePoisson(bk, uk, xkT, xkTdt, nvlong, idefk, ideak, idok);
					}
				      else if(typelong[k] == 3)
					{
					  //add...
					}
				      
				      free(xkTdt);
				  
				      if(asso[5 * jasso + 3] == 1) //linear effect
					{
					  
					  alphakasso = be[nrisq + nes + sumasso] * ytdt;
					  
					  sumasso += 1;
					}
				      else if(asso[5 * jasso + 3] == 2) //polynomial nonlinear effect
					{
					  alphakasso = polynomials(ytdt, be[nrisq + nes + sumasso], be[nrisq + nes + sumasso + 1], be[nrisq + nes + sumasso + 2]);
				      
					  sumasso += 3;
					}
				      else if(asso[5 * jasso + 3] == 3) //logitic nonlinear effect
					{
					  alphakasso = logistic(ytdt, be[nrisq + nes + sumasso], be[nrisq + nes + sumasso + 1], be[nrisq + nes + sumasso + 2]);
				      
					  sumasso += 3;
					}
				    }
			     
				  // interaction with covariate
				  if(asso[5 * jasso + 4] > 0)
				    {
				      alphaAsso += alphakasso * xe[asso[5 * jasso + 4] - 1];
				    }
				  else
				    {
				      alphaAsso += alphakasso;
				    }
				  
				  jasso += 1;
				}
      /* if(l < 2){ */
      /* 	printf("alphaytdt = %f \n ", alphaytdt); */
      /* printf("\n "); */
      /* } */
			      sumnpm2 += npmlong[k];
			      sumnRE += nRE[k];
			      
			      free(bk);
			      free(uk);
			    } // end loop k
			  free(idefk);
			  free(ideak);
			  free(idok);
			  free(xkT);
			} // end if nlong >0
		      // lambda0(T) * exp(Xgamma + sum_k g(yk)*alphak)
		      if((p == 0) & (D == 1)) fevt = baserisk * exp(xegamma + alphaAsso);
		  
      /* if(l < 2){ */
      /* 	printf("fevt = %f \n ", fevt); */
      /* printf("\n "); */
      /* } */
		      // sum_p lambda0(p) * exp(Xgamma + sum_k g(yt)*alphak) * wp
		      if((p > 0) & (p < 16)) cumhaz += baserisk * exp(xegamma + alphaAsso) * wGK[p - 1] * T / 2;

      /* if(l < 2){ */
      /* 	printf("cumhaz = %f \n ", cumhaz); */
      /* printf("\n "); */
      /* } */
		      if(p > 15) trunccumhaz += baserisk * exp(xegamma + alphaAsso) * wGK[p - 16] * T0 / 2;
		      
		    } // end loop p
		  
		  // f(T,D | RE_l)
		  vraisl = vraisl * fevt * exp(-cumhaz);	      
		  
		  // S(T0 | RE_l)
		  if(lefttrunc == 1) truncl = truncl * exp(-trunccumhaz);

		  sumnTi++;
		} // end loop ke
	      
	      sumnpm1 += npmsurv[e];
	      free(be);
	      jasso2 = jasso;
	      
	      if(hazard[e] == 2)
		{
		  sumhazardnodes += nhazardnodes[e];
		  
		  free(nodes);
		}
	    } // end loop e

	  free(idsurve);
	  free(xe);
	} // end survival

      //if((l > 20) & (l < 30)) printf("apres surv vraisl = %f \n", vraisl);
      // loop over the longitudinal models f(Y | RE_l)
      if(nlong > 0)
	{
	  // allocation of indicators of fixed and random effects for marker k
	  idefk = (int *) malloc(nvlong * sizeof(int));
	  ideak = (int *) malloc(nvlong * sizeof(int));
	  idok = (int *) malloc(nvlong * sizeof(int));
  
	  sumni = 0;
	  sumnRE = 0;
	  sumnpm1 = avtlong;
	  for(k = 0; k < nlong; k++)
	    {
	      // observations of marker k
	      if(typelong[k] == 1)
		{
		  yk = (double *) malloc(ni[k] * sizeof(double));
		  for(j = 0; j < ni[k]; j++)
		    yk[j] = Y[sumni + j];
		}
	      else
		{
		  intyk = (int *) malloc(ni[k] * sizeof(int));
		  for(j = 0; j < ni[k]; j++)
		    intyk[j] = (int) (Y[sumni + j] + 0.1);
		}

	      // covariates for marker k
	      xk = (double *) malloc(ni[k] * nvlong * sizeof(double));
	      for(j = 0; j < ni[k]; j++)
		for(m = 0; m < nvlong; m++) // X ligne par ligne
		  {
		    xk[nvlong * j + m] = Xlong[nvlong * (sumni + j) + m];
		  }

	      // parameters for marker k
	      bk = (double *) malloc(npmlong[k] * sizeof(double));
	      for(i = 0; i < npmlong[k]; i++)
		bk[i] = btot[sumnpm1 + i];
	      
	      // indicator of fixed and random effects for marker k
	      for(j = 0; j < nvlong; j++)
		{
		  idefk[j] = idef[k * nvlong + j];
		  ideak[j] = idea[k * nvlong + j];
		}
	      
	      // RE for marker k
	      uk = (double *) malloc(nRE[k] * sizeof(double));
	      for(i = 0; i < nRE[k]; i++)
		uk[i] = U[sumnRE + i];

	      /* if(l < 2){ */
	      /* for(i = 0; i < ni[k]; i++) */
	      /* 	printf("yk= %f \n", yk[i]); */
	      /* for(i = 0; i < npmlong[k]; i++) */
	      /* 	printf("bk= %f \n", bk[i]); */
	      /* for(i = 0; i < nvlong*ni[k]; i++) */
	      /* 	printf("xk= %f \n", xk[i]); */
	      /* for(i = 0; i < nvlong; i++) */
	      /* 	printf("idk= %d %d \n", idefk[i], ideak[i]); */
	      /* for(i = 0; i < nRE[k]; i++) 
	       	printf("uk= %f \n", uk[i]); 
	       	} */
	      // density of marker k
	      vraisk = 0;
	      if(typelong[k] == 1)
		{
		  // Gaussian outcome
		  vraisk = densityGaussian(bk, yk, xk, ni[k], nvlong, idefk, ideak, uk);
		}
	      else if(typelong[k] == 2)
		{
		  // Poisson
		  vraisk = densityPoisson(bk, intyk, xk, ni[k], nvlong, idefk, ideak, idok, uk);
		}
	      else if(typelong[k] == 3)
		{
		  //add...
		}

      /* if(l < 2){ */
      /* 	printf("vraisk = %f \n ", vraisk); */
      /* printf("\n "); */
      /* } */
	      // f(Y | RE_l)
	      vraisl = vraisl * vraisk;
	      
	      sumni += ni[k];
	      sumnpm1 += npmlong[k];
	      sumnRE += nRE[k];

	      if(typelong[k] == 1) free(yk); else free(intyk);
	      free(xk);
	      free(bk);
	      free(uk);
	    } //  end loop k

	  free(idefk);
	  free(ideak);
	  free(idok);
	} // end longitudinal

      //if((l > 20) & (l < 30)) printf("apres long vraisl = %f \n", vraisl);
      *loglik += vraisl / nQMC;
      if(lefttrunc == 1) trunc += truncl / nQMC;
      
      //if((l > 20) & (l < 30)) printf("fin l=%d loglik = %f \n", l, *loglik);

      sumQMC += nREtot;
    } // end QMC

  
  //printf("end before log : %f, log=%f \n ", *loglik, log(*loglik));
  //  printf("\n ");
  // log-likelihood
  if(lefttrunc == 1)
    *loglik = log(*loglik) - log(trunc);
  else
    *loglik = log(*loglik);
  
  // printf("end loglik = %f , lefttrunc=%d \n ", *loglik, lefttrunc);
  //  printf("\n ");
  free(btot);
  free(U);
  free(bRE);
  if(nsurv > 0)
    {
      free(ptsGK);
      free(wGK);
    }

  return;
}
