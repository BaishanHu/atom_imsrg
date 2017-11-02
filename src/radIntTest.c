//#include <boost/math/special_functions/gamma.hpp>
#include <gsl/gsl_sf_gamma.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int
phase(int x)
{
	return x%2==0 ? 1 : -1;
}

double
TalmiI(int p, double k)
{
	  return gsl_sf_gamma(p+1.5+0.5*k) / gsl_sf_gamma(p+1.5);
}

double
TalmiB(int na, int la, int nb, int lb, int p)
{
	if ( (la+lb)%2>0 ) return 0;

	int q = (la+lb)/2;

	double lna = (2*p+1);
	double lnb = (p);
	double lnc = (na);	
	double lnd = (nb);
	double lne = (2*na+2*la+1);
	double lnf = (2*nb+2*lb+1);
	double lng = (na+la);
	double lnh = (nb+lb);
	//double lMax = max(lna, max(lnb, max(lnc, max(lnd, max(lne, max(lnf, max(lng, lnh)))))));

	//boost::multiprecision::float128 B1 = AngMom::phase(p-q);
	double B1 = phase(p-q);
	B1 /= pow(2,(na+nb));

	//B1 *= boost::math::tgamma_ratio( lna+1, lnb+1 );
	//B1 *= sqrt( boost::math::tgamma_ratio( lne+1, lng+1 ) );
	//B1 *= sqrt( boost::math::tgamma_ratio( lnf+1, lnh+1 ) );
	B1 *= gsl_sf_fact(lna) / gsl_sf_fact(lnb);
	B1 *= sqrt( gsl_sf_fact(lne) / gsl_sf_fact(lng) );
	B1 *= sqrt( gsl_sf_fact(lnf) / gsl_sf_fact(lnh) );
	B1 *= sqrt( gsl_sf_fact(lnc) );
	B1 *= sqrt( gsl_sf_fact(lnd) );

	//boost::multiprecision::float128 B2 = 0;
	double B2 = 0;
	int kmin = fmax(0, p-q-nb);
	int kmax = fmin(na, p-q);
	for (int k=kmin;k<=kmax;++k)
	{
		double temp = 1;
		//temp *= ms.GetFactorial(la+k);
		//temp *= ms.GetFactorial(p-int((la-lb)/2)-k);
		//temp /= ms.GetFactorial(k);
		temp /= gsl_sf_gamma(k+1);
		//temp /= ms.GetFactorial(2*la+2*k+1);
		//temp /= ms.GetFactorial(na-k);
		temp /= gsl_sf_gamma(na-k+1);
		//temp /= ms.GetFactorial(2*p-la+lb-2*k+1);
		//temp /= ms.GetFactorial(nb - p + q + k);
		temp /= gsl_sf_gamma(nb - p + q + k + 1);
		//temp /= ms.GetFactorial(p-q-k);
		temp /= gsl_sf_gamma(p-q-k+1);
		//temp *= boost::math::tgamma_ratio( la+k +1, 2*la+2*k+1 +1);
		temp *= gsl_sf_fact( la + k ) / gsl_sf_fact( 2*la + 2*k + 1 );
		temp *= gsl_sf_fact( p-((la-lb)/2)-k ) / gsl_sf_fact( 2*p-la+lb-2*k+1 );
		//temp *= boost::math::tgamma_ratio( p-int((la-lb)/2)-k +1, 2*p-la+lb-2*k+1 +1);
		B2 += temp; 
	}
	double retB = (B1 * B2);//.convert_to<double>();
	return retB; //B1 * B2;
}

double
RadialIntegral_RpowK(int na, int la, int nb, int lb, int k)
{
	double I = 0;
	/* mpf_t I, tb, ti;
	mpf_init2 (I, 256);
	mpf_init2 (tb, 256);
	mpf_init2 (ti, 256);
	I = 0; */
	int pmin = (la+lb)/2;
	int pmax = pmin + na + nb;
	for (int p=pmin;p<=pmax;++p)
	{
		double tb = TalmiB(na,la,nb,lb,p);
		double ti = TalmiI(p,k);
		//printf ("p      = %d\n", p);
		//printf ("TalmiB = %.18f\n", tb);
		//printf ("TalmiI = %.18f\n", ti);
		I += tb * ti;
	}
	//mpf_clear (tb);
	//mpf_clear (ti);
	return I;
} 

double
RadialIntegral(int na, int la, int nb, int lb, int L)
{
if ((la+lb+L)%2!=0) return RadialIntegral_RpowK(na,la,nb,lb,L);
int tau_a = fmax((lb-la+L)/2,0);
int tau_b = fmax((la-lb+L)/2,0);
int sigma_min = fmax(fmax(na-tau_a,nb-tau_b),0);
int sigma_max = fmin(na,nb);

double term1 = phase(na+nb) * gsl_sf_fact(tau_a)*gsl_sf_fact(tau_b) * sqrt(gsl_sf_fact(na)*gsl_sf_fact(nb)
                   / (gsl_sf_gamma(na+la+1.5)*gsl_sf_gamma(nb+lb+1.5) ) );
double term2 = 0;
for (int sigma=sigma_min; sigma<=sigma_max; ++sigma)
{
term2 += gsl_sf_gamma(0.5*(la+lb+L)+sigma+1.5) / (gsl_sf_fact(sigma)*gsl_sf_fact(na-sigma)*gsl_sf_fact(nb-sigma)*gsl_sf_fact(sigma+tau_a-na)*gsl_sf_fact(sigma+tau_b-nb) );
}
return term1*term2;

}

int
main()
{
	//int na = 24;
	int la = 0;
	//int nb = 24;
	int lb = 0;
	int L = -1;
	int nstart = 20;
	for ( int na = nstart; na <= 25; na++)
	{
		for ( int nb = nstart; nb <= na; nb ++)
		{
			printf ( "+++++++++++++++++++++++++++++++\n" );
			printf ( "na      = % d\n", na);
			printf ( "nb      = % d\n", nb);
			printf ( "la      = % d\n", la);
			printf ( "lb      = % d\n", lb);
			printf ( "L       = % d\n", L);
			double rad = RadialIntegral( na,la, nb,lb, L);
			//cout << "Rad =" << rad << endl;
			//if ( abs(rad) > 1.1
			printf ( "Rad     = % .18f\n", rad );
			printf ( "\n" );
		}
	}

	return 0;
}
