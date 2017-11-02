/*
   Example adapted from the GNU Scientific Library Reference Manual
    Edition 1.1, for GSL Version 1.1
    9 January 2002
   URL: gsl/ref/gsl-ref_23.html#SEC364
    Revised: 19-Nov-2002 by Dick Furnstahl  furnstahl.1@osu.edu
*/

/* 
  Compile and link with:
    gcc -c qagiu_test.c
    gcc -o qagiu_test qagiu_test.o -lgsl -lgslcblas -lm
*/

/* 
   Each algorithm computes an approximation to a definite integral of
   the form,

   I = \int_a^b f(x) w(x) dx

   where w(x) is a weight function (for general integrands w(x)=1). The
   user provides absolute and relative error bounds (epsabs, epsrel)
   which specify the following accuracy requirement,

   |RESULT - I|  <= max(epsabs, epsrel |I|)

   where RESULT is the numerical approximation obtained by the
   algorithm. The algorithms attempt to estimate the absolute error
   ABSERR = |RESULT - I| in such a way that the following inequality
   holds,

   |RESULT - I| <= ABSERR <= max(epsabs, epsrel |I|)

   The routines will fail to converge if the error bounds are too
   stringent, but always return the best approximation obtained up to
   that stage.
*/

/*
   QAGIU adaptive integration on infinite intervals

   Function: 
   int gsl_integration_qagiu (gsl_function * f, double a,
      double epsabs, double epsrel, size_t limit, 
      gsl_integration_workspace * workspace, 
      double *result, double *abserr)

   This function computes the integral of the function f over the
   semi-infinite interval (a,+\infty). The integral is mapped onto the
   interval (0,1] using the transformation x = a + (1-t)/t,

   \int_{a}^{+\infty} dx f(x) = 
        \int_0^1 dt f(a + (1-t)/t)/t^2

   and then integrated using the QAGS algorithm.
*/


/*
   The integrator QAGIU will handle a large class of infinite-range integrals.
   For example, consider the following integral, 
   
   \int_1^\infty cos[2x] e^{-x} dx = -0.16442310483055015762

   The program below computes this integral to a relative accuracy bound
   of 1e-8.
*/

/*
   The results below show that the desired accuracy is achieved after
   10 subdivisions.

   result          = -0.164423104830549505
   exact result    = -0.164423104830550171
   estimated error =  0.000000005882847412
   actual error    =  0.000000000000000666
   intervals =  10
*/

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_laguerre.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>

double
my_integrand (double x, void *params_ptr)
{
  double alpha;

  alpha = *(double *) params_ptr;

  const int n = 1;
  printf ("n          = % .18f\n", n);
  return (cos (alpha * n * x) * exp (-x));
}

struct my_f_params {int n; int l; int np; int Z;};

double
HO( double x, void * p )
{
	struct my_f_params * params = (struct my_f_params *)p;
        int n = (params->n);
        int l = (params->l);
	int np = (params->np);
	int Z = (params->Z);
	double a = 0.0529; // Bohr Radius 0.53 A, 0.0529 nm
	double ho = 0;
	double m = 1;
	double h = 1;
	double w = 2*13.605;//13.605 * 2 / h; // wavelength of oscillator 13.605 *2 ? 
	double v = m*w/(2*h);
	int test = 0;
	double c = Z / (n * a);
	double cp = Z / (np * a);

	double hy = sqrt( pow(2*c,3) * gsl_sf_fact(n-l-1)/((2*n*gsl_sf_fact(n+l)) ) ) * pow(2*c,l) * exp( -c*x ) * gsl_sf_laguerre_n(n-l-1, 2*l+1, 2*x*c);
	//double hyp = sqrt( pow(2*cp,3) * gsl_sf_fact(np-l-1)/((2*n*gsl_sf_fact(np+l)) ) ) * pow(2*cp,l) * exp( -cp*x ) * gsl_sf_laguerre_n(np-l-1, 2*l+1, 2*x*cp);
	//double o = sqrt( sqrt( 2* pow(13.605,3) / 3.14159 ) * pow(2, n+2*l+3) * gsl_sf_fact(n) * pow(13.605,l) / gsl_sf_doublefact(2*n+2*l+1) ) * exp( -v*x*x ) * gsl_sf_laguerre_n(n, l+0.5, 2*v*x*x);
	double op = sqrt( sqrt( 2* pow(13.605,3) / 3.14159 ) * pow(2, np+2*l+3) * gsl_sf_fact(np) * pow(13.605,l) / gsl_sf_doublefact(2*np+2*l+1) ) * exp( -v*x*x ) * gsl_sf_laguerre_n(np, l+0.5, 2*v*x*x);
	ho = hy*op;
	//ho = hy*hyp; // Correctly orthonormal!
	//ho = o*op; // Correctly orthonormal! 
	//ho = o*o;
	//ho = hy*hy;
	//int l = 0;
	//ho = no * nhy * pow( x, 2*l+2 ) * exp( -v*x*x - x/(ny * a) ) * gsl_sf_laguerre_n(n, l+0.5, 2*v*x*x) * gsl_sf_laguerre_n(ny-l-1, 2*l+1, 2*x/(ny * a));
	return ho*x*x; // x*x comes from integrating of r^2 dr
}



double
OsToHydroCoeff( double x, void * p )
{
	struct my_f_params * params = (struct my_f_params *)p;
        int n = (params->n);
        int l = (params->l);
	int np = (params->np);
	int Z = (params->Z);
	double c = Z / (n * 0.0529); // 
	double m = 1; 			// Electron mass, in atomic units
	double h = 1; 			// Reduced plancks' constant, in atomic units
	double w = 13.605 * 2 / h; 	// wavelength of oscillator 13.605 *2 ? 
	double v = m*w/(2*h);
	       
	return pow( x, 2*l+2 ) * exp( -v*x*x - c*x ) * gsl_sf_laguerre_n(np, l+0.5, 2*v*x*x) * gsl_sf_laguerre_n(n-l-1, 2*l+1, 2*x*c);
}

int
main (void)
{
  int size = 1000;
  gsl_integration_workspace *work_ptr =
    gsl_integration_workspace_alloc (size);

  gsl_set_error_handler_off ();

  double lower_limit = 0;		/* start integral from lower_limit (to infinity) 	*/
  long double abs_error = 1.0e-9;	/* to avoid round-off problems 				*/
  long double rel_error = 1.0e-9;	/* the result will usually be much better 		*/
  double result=0;			/* the result from the integration 			*/
  double error;				/* the estimated error from the integration 		*/
  double temp = 0;
  double eH = -4/pow(2,2);
  double tmax = 0;
  double prev = 0;
  struct my_f_params alpha = {1,0,1,1};
  gsl_function My_function;
  void *params_ptr = &alpha;
  int Z = 2;

  for (int n = 1; n < 2; n++)
  {
   for (int l = 0; l < 1; l++)
   {
    double hydrogenCoeff = sqrt( pow(2*Z/(n * 0.0529),3) * gsl_sf_fact(n-l-1)/((2*n*gsl_sf_fact(n+l)) ) ) * pow(2*Z/(n * 0.0529),l);
    for (int np = 0; np < 100; np++)
    {
     double OscilCoeff = sqrt( sqrt( 2* pow(13.605,3) / 3.14159 ) * pow(2, np+2*l+3) * gsl_sf_fact(np) * pow(13.605,l) / gsl_sf_doublefact(2*np+2*l+1) );
     result = 0;
     //printf ("n                 = % .18f\n", n);
     //printf ("l                 = % .18f\n", l);
     //printf ("np                = % .18f\n", np);
     //double alpha = my_f_params(0,0,0);
     alpha.n = n;
     //printf ("alpha.n         = % d\n", alpha.n);
     alpha.l = l;
     //printf ("alpha.l         = % d\n", alpha.l);
     alpha.np = np;
     //printf ("alpha.np        = % d\n", alpha.np);
     alpha.Z = Z;

     //My_function.function = &OsToHydroCoeff;
     My_function.function = &HO;
     //My_function.function = &my_integrand;
     My_function.params = params_ptr;
     //My_function.function->n = n;

     gsl_integration_qagiu (&My_function,
			lower_limit,
			abs_error,
			rel_error,
			size,
			work_ptr,
			&result,
			&error);
     //temp += result*result*(2*alpha.np + alpha.l + 3./2)*2;
     //if ( sqrt(pow(temp,2)) > tmax ) tmax = temp;
     //printf ("tmax            = % .18f\n", tmax);

     //printf ("percent diff    = % .18f\n", sqrt(pow(temp-tmax,2))/tmax);
     //if ( (temp <= 1.03*eH) && (temp >= 0.97*eH) ) printf ("Here is one.\n");
     //printf ("Energy          = % .18f\n", temp);
     //printf ("Prev Energy     = % .18f\n", prev);
     //prev = temp;
     //const double alpha[3]={0.1, 0.1, 1}, theta[3]={1/3., 1/3., 1/3.};
     //double x;
     //x=gsl_ran_dirichlet_pdf(3, alpha, theta);
     //printf ("x               = % .18f\n", x);
     //printf ("ln(gamma(0))    = % .18f\n", gsl_sf_lngamma(0));
     //printf ("ln(gamma(1))    = % .18f\n", gsl_sf_lngamma(1));
     printf ("%.18f\n", result);
     //printf ("%.18f\n", result*result);
     result = 0;
     //printf ("estimated error = % .18f\n", error);
     //printf ("intervals       =  %d\n\n", work_ptr->size);
    }
   }
  }
  return 0;
}
