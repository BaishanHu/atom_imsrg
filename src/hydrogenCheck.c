#include <stdio.h>
#include <gsl/gsl_sf_laguerre.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_math.h>
#include <math.h>

#include "hydrogenCheck.h"

/*
double H( double x, void * params )
{
	double a = 1; // Bohr Radius
	int n = *(int *) params;
	//long double h = sqrt( pow(2/(n*a),3) ) * gsl_fact(n-l-1)/gsl_fact(2*n*(n+l)) ) * exp(-r/(n*a)) * pow(2*r/(n*a)),l) * gsl_sf_laguerre_n(n, 2*l+1, x);
	long double h = sqrt( pow(2/(n*a),3) ) * gsl_sf_fact(n-1)/( 2*n*gsl_sf_fact((n)) ) * exp(-x/(n*a)) * gsl_sf_laguerre_n(n, 1, x);
	return h;
}

double O( double x, void * params )
{
	int n = *(int *) params;
	double m = 1; // mass of the oscillator, should be electron mass
	double w = 27.4113;
	double h = 1; // hbar - reduced Planck constant
	double v = m*w/(2*h);
	//long double o = sqrt( sqrt( 2*(v)^3/pi * 2^(n+2*l+3)*gsl_fact(n)*pow(v,l) / gsl_doublefact(2*n+2*l+1) ) * pow(x,l) * exp( -v*x^2 ) * gsl_sf_laguerre_n(n, l+0.5, 2*v*x*x);
	long double o = sqrt( sqrt( 2*pow(v,3)/3.14159 * pow(2,(n+3))*gsl_sf_fact(n) / gsl_sf_doublefact(2*n+1) ) * exp( -v*x*x ) * gsl_sf_laguerre_n(n, 0.5, 2*v*x*x));
	return o;
}  */

double HO( double x, void * params )
{
	double a = 1; // Bohr Radius
	int n = *(int *) params;
	//long double h = sqrt( pow(2/(n*a),3) ) * gsl_fact(n-l-1)/gsl_fact(2*n*(n+l)) ) * exp(-r/(n*a)) * pow(2*r/(n*a)),l) * gsl_sf_laguerre_n(n, 2*l+1, x);
	long double hy = sqrt( pow(2/(1*a),3) ) * gsl_sf_fact(1-1)/( 2*1*gsl_sf_fact((1)) ) * exp(-x/(1*a)) * gsl_sf_laguerre_n(1, 1, x);
	double m = 1; // mass of the oscillator, should be electron mass
	double w = 27.4113;
	double h = 1; // hbar - reduced Planck constant
	double v = m*w/(2*h);
	//long double o = sqrt( sqrt( 2*(v)^3/pi * 2^(n+2*l+3)*gsl_fact(n)*pow(v,l) / gsl_doublefact(2*n+2*l+1) ) * pow(x,l) * exp( -v*x^2 ) * gsl_sf_laguerre_n(n, l+0.5, 2*v*x*x);
	long double o = sqrt( sqrt( 2*pow(v,3)/3.14159 * pow(2,(n+3))*gsl_sf_fact(n) / gsl_sf_doublefact(2*n+1) ) * exp( -v*x*x ) * gsl_sf_laguerre_n(n, 0.5, 2*v*x*x));
	double ho = hy*o;
	return ho;
}

int main (void)
{
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

	double lower_limit = 0.;
	double abs_error = 1.0e-8;
	double rel_error = 1.0e-8;
	double result;
	double error;
	double n=1;

	gsl_function ho;
	void *params = &n;

	ho.function = &HO;
	ho.params = params;

	//for (int n=0; n<5; n++)
	//{
		 //(gsl_function * f, double a, double epsabs, double epsrel, size_t limit, gsl_integration_workspace * workspace, double * result, double * abserr)
		gsl_integration_qagiu (&ho,
					lower_limit,
					abs_error,
					rel_error,
					1000,
					w,
					&result,
					&error);

		printf ("n=		% .18f\n", n);
		printf ("Result=	% .18f\n", result);
		printf ("Error=		% .18f\n", error);
		//gsl_integration_qags (&f, 0, inf
	//}

	gsl_integration_workspace_free (w);

	return 0;
} 
