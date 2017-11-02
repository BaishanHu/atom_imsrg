#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_laguerre.h>
#include <gsl/gsl_math.h>

struct my_f_params {int n; int l; int np;};

double
OsToHydroCoeff( double x, void * p )
{
	struct my_f_params * params = (struct my_f_params *)p;
        int n = (params->n);
        int l = (params->l);
	int np = (params->np);
	double c = 1 / (n * 0.0529);
	double w = 13.605;
	       
	return pow( x, 2*l+2 ) * exp( -w*x*x - c*x ) * gsl_sf_laguerre_n(np, l+0.5, 2*w*x*x) * gsl_sf_laguerre_n(n-l-1, 2*l+1, 2*w*c);
}

int
main (void)
{
  int size = 1000;
  gsl_integration_workspace *work_ptr =
    gsl_integration_workspace_alloc (size);

  gsl_set_error_handler_off ();

  double lower_limit = 0;		/* start integral from lower_limit (to infinity) 	*/
  long double abs_error = 1.0e-6;	/* to avoid round-off problems 				*/
  long double rel_error = 1.0e-6;	/* the result will usually be much better 		*/
  double result=0;			/* the result from the integration 			*/
  double error;				/* the estimated error from the integration 		*/

  struct my_f_params alpha = {1,0,1};
  gsl_function My_function;
  void *params_ptr = &alpha;

  for (int n = 1; n < 6; n++)
  {
   for (int l = 0; l < n; l++)
   {
    for (int np = 0; np < 6; np++)
    {
     result = 0;
     alpha.n = n;
     printf ("alpha.n         = % d\n", alpha.n);
     alpha.l = l;
     printf ("alpha.l         = % d\n", alpha.l);
     alpha.np = np;
     printf ("alpha.np        = % d\n", alpha.np);


     My_function.function = &OsToHydroCoeff;
     My_function.params = params_ptr;

     gsl_integration_qagiu (&My_function,
			lower_limit,
			abs_error,
			rel_error,
			size,
			work_ptr,
			&result,
			&error);

     printf ("result         =% .18f\n", result);
     result = 0;

    }
   }
  }
  return 0;
}
