#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include <stdio.h>
#include <gsl/gsl_sf_laguerre.h>
#include <gsl/gsl_sf_gamma.h>

/* Computation of the integral,

      I = int (dx dy dz)/(2pi)^3  1/(1-cos(x)cos(y)cos(z))

   over (-pi,-pi,-pi) to (+pi, +pi, +pi).  The exact answer
   is Gamma(1/4)^4/(4 pi^3).  This example is taken from
   C.Itzykson, J.M.Drouffe, "Statistical Field Theory -
   Volume 1", Section 1.1, p21, which cites the original
   paper M.L.Glasser, I.J.Zucker, Proc.Natl.Acad.Sci.USA 74
   1800 (1977) */

/* For simplicity we compute the integral over the region 
   (0,0,0) -> (pi,pi,pi) and multiply by 8 */

double exact = 27.2114;//0.0529*3/2;

struct my_f_params { int n1; int l1; int n2; int l2; int n3; int l3; int n4; int l4; int Z;};

int lim = 1;
double a = 0.0529; // Bohr Radius 0.53 A, 0.0529 nm
int HBARC = 197.3; // 197.3 eV nm or 1973 ev A

/*
Integral of <12,J|1/abs(r3-r4)|34,J>
*/
double
e2b (double *k, size_t dim, void * p)
{
  (void)(dim); /* avoid unused parameter warnings */
  struct my_f_params * params = (struct my_f_params *)p;

  double x3 = k[0]/(lim-k[0]);
  double x4 = k[1]/(lim-k[1]);

  int Z = (params->Z);

  int n1 = (params->n1);
  int l1 = (params->l1);
  double c1 = Z / (n1 * a);
  double h1 = pow(2*c1,3) * gsl_sf_fact(n1-l1-1) / ( 2*n1*gsl_sf_fact(n1+l1) );
  double hy1 = sqrt( h1 ) * pow(2*c1*x3,l1) * exp( -c1*x3 ) * gsl_sf_laguerre_n(n1-l1-1, 2*l1+1, 2*x3*c1);

  int n2 = (params->n2);
  int l2 = (params->l2);
  double c2 = Z / (n2 * a);
  double h2 = pow(2*c2,3) * gsl_sf_fact(n2-l2-1) / ( 2*n2*gsl_sf_fact(n2+l2) );
  double hy2 = sqrt( h2 ) * pow(2*c2*x4,l2) * exp( -c2*x4 ) * gsl_sf_laguerre_n(n2-l2-1, 2*l2+1, 2*x4*c2);

  int n3 = (params->n3);
  int l3 = (params->l3);
  double c3 = Z / (n3 * a);
  double h3 = pow(2*c3,3) * gsl_sf_fact(n3-l3-1) / ( 2*n3*gsl_sf_fact(n3+l3) );
  double hy3 = sqrt( h3 ) * pow(2*c3*x3,l3) * exp( -c3*x3 ) * gsl_sf_laguerre_n(n3-l3-1, 2*l3+1, 2*x3*c3);

  int n4 = (params->n4);
  int l4 = (params->l4);
  double c4 = Z / (n4 * a);
  double h4 = pow(2*c4,3) * gsl_sf_fact(n4-l4-1) / ( 2*n4*gsl_sf_fact(n4+l4) );
  double hy4 = sqrt( h4 ) * pow(2*c4*x4,l4) * exp( -c4*x4 ) * gsl_sf_laguerre_n(n4-l4-1, 2*l4+1, 2*x4*c4);
  //double den = 

  double A = hy1 * hy2 * 1/sqrt( fabs( pow(x4,2) + pow(x3,2) -2*x3*x4*cos(x3) ) ) * hy3 * hy4 * 1/pow(lim-k[0],2) * 1/pow(lim-k[1],2) * x3*x3 * x4*x4 * HBARC/137; // hbarc*alpha
  return A;
}

double
eb (double *k, size_t dim, struct my_f_params * p ) //void * p)
{
  (void)(dim); /* avoid unused parameter warnings */
  struct my_f_params * params = (struct my_f_params *)p;

  double x3 = k[0]/(lim-k[0]);
  double x4 = k[1]/(lim-k[1]);

  int Z = (params->Z);

  int n1 = (params->n1);
  int l1 = (params->l1);
  double c1 = Z / (n1 * a);
  double h1 = pow(2*c1,3) * gsl_sf_fact(n1-l1-1) / ( 2*n1*gsl_sf_fact(n1+l1) );
  double hy1 = sqrt( h1 ) * pow(2*c1*x3,l1) * exp( -c1*x3 ) * gsl_sf_laguerre_n(n1-l1-1, 2*l1+1, 2*x3*c1);

  int n2 = (params->n2);
  int l2 = (params->l2);
  double c2 = Z / (n2 * a);
  double h2 = pow(2*c2,3) * gsl_sf_fact(n2-l2-1) / ( 2*n2*gsl_sf_fact(n2+l2) );
  double hy2 = sqrt( h2 ) * pow(2*c2*x4,l2) * exp( -c2*x4 ) * gsl_sf_laguerre_n(n2-l2-1, 2*l2+1, 2*x4*c2);

  int n3 = (params->n3);
  int l3 = (params->l3);
  double c3 = Z / (n3 * a);
  double h3 = pow(2*c3,3) * gsl_sf_fact(n3-l3-1) / ( 2*n3*gsl_sf_fact(n3+l3) );
  double hy3 = sqrt( h3 ) * pow(2*c3*x3,l3) * exp( -c3*x3 ) * gsl_sf_laguerre_n(n3-l3-1, 2*l3+1, 2*x3*c3);

  int n4 = (params->n4);
  int l4 = (params->l4);
  double c4 = Z / (n4 * a);
  double h4 = pow(2*c4,3) * gsl_sf_fact(n4-l4-1) / ( 2*n4*gsl_sf_fact(n4+l4) );
  double hy4 = sqrt( h4 ) * pow(2*c4*x4,l4) * exp( -c4*x4 ) * gsl_sf_laguerre_n(n4-l4-1, 2*l4+1, 2*x4*c4);

  double A = hy1 * hy2 * hy3 * hy4 * 1/pow(lim-k[1],2) * 1/pow(lim-k[0],2) * x3*x3 * x4*x4; // hbarc*alpha
  return A;
}

double
g (double *k, size_t dim, void *params)
{
  (void)(dim); /* avoid unused parameter warnings */
  (void)(params);
  double A = 1.0 / (M_PI * M_PI * M_PI);
  return A / (1.0 - cos (k[0]) * cos (k[1]) * cos (k[2]));
}

void
display_results (char *title, double result, double error)
{
  printf ("%s ==================\n", title);
  printf ("result = % .6f\n", result);
  printf ("sigma  = % .6f\n", error);
  printf ("exact  = % .6f\n", exact);
  printf ("error  = % .6f = %.2g sigma\n", result - exact,
          fabs (result - exact) / error);
}

double
integrate ( double N, double M, void * p )
{
  printf ("Entering integrate.\n");
  struct my_f_params * params = (struct my_f_params *)p;
/*
  int Z = (params->Z);

  int n1 = (params->n1);
  int l1 = (params->l1);
  double c1 = Z / (n1 * a);

  int n2 = (params->n2);
  int l2 = (params->l2);
  double c2 = Z / (n2 * a);

  int n3 = (params->n3);
  int l3 = (params->l3);
  double c3 = Z / (n3 * a);

  int n4 = (params->n4);
  int l4 = (params->l4);
  double c4 = Z / (n4 * a); */
  int Z = 2;
  struct my_f_params alpha = { 2,0, 1,0, 1,0, 1,0, Z };
  double dx = 1 / N;
  double dy = 1 / M;
  double k[2] = {0, 0};
  double I = 0;
  for ( double xi = 0; xi < 1; xi += dx )
  {
    for ( double yj = 0; yj < 1; yj += dy )
    {
      //printf ("xi = %.18f\n", xi);
      //printf ("yj = %.18f\n", xi);
      k[0] = xi;
      k[1] = yj;
      I += eb ( k, 0, &alpha );
      //printf ("I = %.18f\n", xi);
    }
  }
  printf ("leaving integrate.\n");
  return I;
}

int
main (void)
{
  double res, err;

  double xl[2] = { 0, 0 }; // Lower limit 0; [0,1)
  double xu[2] = { lim, lim } ; // Upper limit 1;

  const gsl_rng_type *T;
  gsl_rng *r;
  int Z = 2;
  struct my_f_params alpha = { 1,0, 1,0, 1,0, 1,0, Z };
  gsl_monte_function G = { &eb, 2, &alpha };

  size_t calls = 1e6;

  gsl_rng_env_setup ();

  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
  
  double self = integrate ( 1000, 1000, &alpha );
  printf ("self = %.18f\n", self);

  {
    gsl_monte_plain_state *s = gsl_monte_plain_alloc (2);
    gsl_monte_plain_integrate (&G, xl, xu, 2, calls, r, s, 
                               &res, &err);
    gsl_monte_plain_free (s);

    display_results ("plain", res, err);
  }

  {
    gsl_monte_miser_state *s = gsl_monte_miser_alloc (2);
    gsl_monte_miser_integrate (&G, xl, xu, 2, calls, r, s,
                               &res, &err);
    gsl_monte_miser_free (s);

    display_results ("miser", res, err);
  }

  {
    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (2);
    gsl_monte_vegas_integrate (&G, xl, xu, 2, calls/5, r, s,
                               &res, &err);
    display_results ("vegas warm-up", res, err);
    printf ("converging...\n");

    do
      {
        gsl_monte_vegas_integrate (&G, xl, xu, 2, calls/5, r, s,
                                   &res, &err);
        printf ("result = % .6f sigma = % .6f "
                "chisq/dof = %.1f\n", res, err, gsl_monte_vegas_chisq (s));
      }
    while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.1);

    display_results ("vegas final", res, err);

    gsl_monte_vegas_free (s);
  }

  gsl_rng_free (r);

  return 0;
}
