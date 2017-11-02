#include <stdio.h>
#include <gsl/gsl_sf_laguerre.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_gamma.h>
#include <math.h>

//double H( double x, void * params );
//double O( double x, void * params );
double HO( double x, void * params );

