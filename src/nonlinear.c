#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "nonlinear.h"

// degre 3 polynomials for non linear effect
double polynomials(double y, double b1, double b2, double b3)
{
  double res = 0.0;

  res = y * b1 + pow(y, 2) * b2 + pow(y, 3) * b3 ;

  return res;
}


// logistic for non linear effect with S shape
double logistic(double y, double b1, double b2, double b3)
{
  double res = 0.0;

  res = b1 / (1 + exp(-b2 * (y - b3)));
 
  return res;
}
