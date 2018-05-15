/*
#include "cube_test.hh"
#include <string>
#include <iostream>
#include "cubature.h"
#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_laguerre.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_gamma.h>
*/

int main();
double Yml (double theta, int l, int m);
double hydrogenWF(double x, double theta, int n, int l, int m, int Z);
int f(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);
