#include "cube_test.hh"
#include <string>
#include <iostream>
#include "cubature.h"
#include "hcubature.c"
#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_laguerre.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_gamma.h>

using namespace std;

struct my_f_params { int n1; int l1; int m1; int n2; int l2; int m2; int n3; int l3; int m3; 
int n4; int l4; int m4; int Z;};

double Yml (double theta, int l, int m)
{
	double temp=pow(-1,m) * sqrt( (2*l+1) * gsl_sf_fact(l-m) / (4*3.141592 * gsl_sf_fact(l+m)));
	return temp*gsl_sf_legendre_Plm( l, m, cos(theta) );
}

double hydrogenWF(double x, double theta, int n, int l, int m, int Z)
{
	double a = 0.0529;
        double c = Z/(n * a);
        double h = pow(2*c, 3) * gsl_sf_fact(n-l-1) / (2*n*gsl_sf_fact(n+l) );
        double hy = sqrt(h) * pow(2*c*x, l) * exp(-c*x) * gsl_sf_laguerre_n(n-l-1, 2*l+1, 2*x*c) * Yml(theta, l, m);
	
	return hy;
}

double Rnl(double x, int n, int l, int Z)
{
	n = n+1;
	double a = 1; //0.0529;
        double c = Z/(n * a);
        double h = pow(2*c, 3) * gsl_sf_fact(n-l-1) / (2*n*gsl_sf_fact(n+l) );
	double hy = sqrt(h) * pow(2*c*x, l) * exp(-c*x) * gsl_sf_laguerre_n(n-l-1, 2*l+1, 2*x*c);
	return hy;
}

struct cs_RabRcd_params { int na; int la;
                          int nb; int lb;
                          int nc; int lc;
                          int nd; int ld; int lp; double b; };

double cs_Rnl(double r, int n, int l, double b)
{
	// b = 1;
        double norm = 2./b * gsl_sf_fact(n) * 1./(gsl_sf_fact(n+2*l+2)) ;
        return sqrt(norm) * pow(2.*r/b, l+1) * exp(-r/b) * gsl_sf_laguerre_n(n, 2*l+2, 2.*r/b);
}

double r_exp(double r, int n, int l, double b)
{
	return exp(-(r/b));
}

int cs_RabRcd(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{
    (void)(ndim); // Avoid unused parameter warnings;
    (void)(fval);
    struct cs_RabRcd_params * params = (struct cs_RabRcd_params*)fdata;

    int na = (params->na);
    int nb = (params->nb);
    int nc = (params->nc);
    int nd = (params->nd);

    int la = (params->la);
    int lb = (params->lb);
    int lc = (params->lc);
    int ld = (params->ld);

    int lp = (params->lp);

    double b = (params->b);

    double r1 = x[0]; //*1./(1.-x[0]);
    double r2 = x[1]; //*1./(1.-x[1]);

    double rmin = min(r1,r2); // Factor of 1./b accounted for in CSTwoBody.
    double rmax = max(r1,r2);
/*
    double Ra = r1*Rnl(r1,na,la,b);
    double Rb = r2*Rnl(r2,nb,lb,b);
    double Rc = r1*Rnl(r1,nc,lc,b);
    double Rd = r2*Rnl(r2,nd,ld,b);
*/
    double Ra = cs_Rnl(r1,na,la,b);
    double Rb = cs_Rnl(r2,nb,lb,b);
    double Rc = cs_Rnl(r1,nc,lc,b);
    double Rd = cs_Rnl(r2,nd,ld,b);

    double Ha = 27.2114;

    // fval[0] = Ra*Rb*Rc*Rd * pow(rmin,lp) *1./ pow(rmax,lp+1);// * 1./pow(1-x[0],2) * 1./pow(1-x[1],2);

    fval[0] = Ra*Rb*Rc*Rd/(r1);// * 1./pow(1-x[0],2) * 1./pow(1-x[1],2);

    return 0;
}

//struct cs_RabRcd_params int_params = { o1.n,o1.l, o2.n,o2.l, o3.n,o3.l, o4.n,o4.l, L, b};
//double val;
//double err;
//hcubature(1, &cs_RabRcd, &int_params, 2, xmin, xmax, max_iter, 0, max_err, ERROR_INDIVIDUAL, &val, &err);


int f(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval) {
	(void)(ndim); // Avoid unused parameter warnings;
	(void)(fval);
	struct my_f_params * params = (struct my_f_params *)fdata;
	int lim = 1;
	float a = 0.0529; // Bohr Radius
	double x3 = x[0]/(1-x[0]);
	double x4 = x[2]/(1-x[2]);
	double theta3 = x[1];
	double theta4 = x[3];
	int Z = (params->Z);
	double del = 1e-9;
	double PI = 3.141592;

	// First function
	int n1 = (params->n1);
	int l1 = (params->l1);
	int m1 = (params->m1);
	
	// Second Function
        int n2 = (params->n2);
        int l2 = (params->l2);
        int m2 = (params->m2);
 
        // Third function
        int n3 = (params->n3);
        int l3 = (params->l3);
        int m3 = (params->m3);

        // Fourth Function
        int n4 = (params->n4);
        int l4 = (params->l4);
        int m4 = (params->m4);

	// The phi components are orth
	if (m1 != m3 || m2 != m4)
		return 0;
	// l < n
	if (l1 >= n1 || l2 >= n2 || l3 >= n3 || l4 >= n4)
		return 0;

	// -l <= m <= l 
	//if (m1 > l1 || m2 > l2 || m3 > l3 || m4 > l4 )
	//	return 0;

	double h1 = hydrogenWF(x3, theta3, n1, l1, m1, Z);
        double h2 = hydrogenWF(x4, theta4, n2, l2, m2, Z);
        double h3 = hydrogenWF(x3, theta3, n3, l3, m3, Z);
        double h4 = hydrogenWF(x4, theta3, n4, l4, m4, Z);

	fval[0] = h1 * h3 * x3*x3 * 1/pow(1-x[0],2);
	fval[0]*= h2 * h4 * x4*x4 * 1/pow(1-x[2],2);
	fval[0]*= sin(theta3) * sin(theta4) * 4*PI*PI;
	fval[0]*= 1/sqrt( x3*x3 + x4*x4 - 2*x3*x4*cos(theta3) );

	return 0;
}

int main()
{
	double xmin[2] = {0,0};
	double pi = 3.141592;
	double hbarc = 197.32697;
	double fs = 1/137.035999139;
	double xmax[2] = {1000,1000};
	double val = 0;
	double err = 0;
	double b = 1;
	struct cs_RabRcd_params alpha = { 0,0, 0,0, 0,0, 0,0, 0, b};
	hcubature(1, &cs_RabRcd, &alpha, 2, xmin, xmax, 1e6, 0, 1e-5, ERROR_INDIVIDUAL, &val, &err);
	cout << "b=" << b << endl;
	cout << "The result is: " << val << " in eV: " << val*27.2114 << endl;
	cout << "The error is: " << err << endl;
	return 0;
}
