#ifndef AUXFCT_H_   
#define AUXFCT_H_


void matmulC (double * m1, double * m2, int l1, int n1, int c2, double * res);

void matmul (double * m1, double * m2, int l1, int n1, int c2, double * res);

void ltmul(double * mat, double * v, int n);

void linpred(double * beta, double * X, double * u, int * idef, int * idea, int ni, int nv, double * res);

#endif
