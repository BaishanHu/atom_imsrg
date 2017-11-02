
#include "imsrg_util.hh"
#include "AngMom.hh"
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/multiprecision/float128.hpp>
#include <boost/implicit_cast.hpp>
#include <gsl/gsl_integration.h>
#include <list>
#include <cmath>
#include <algorithm>

using namespace AngMom;

namespace imsrg_util
{

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
		double lMax = max(lna, max(lnb, max(lnc, max(lnd, max(lne, max(lnf, max(lng, lnh)))))));

		//boost::multiprecision::float128 B1 = AngMom::phase(p-q);
		long double B1 = phase(p-q);
		B1 /= pow(2,(na+nb));

		B1 *= boost::math::tgamma_ratio( lna+1, lnb+1 );
		B1 *= sqrt( boost::math::tgamma_ratio( lne+1, lng+1 ) );
		B1 *= sqrt( boost::math::tgamma_ratio( lnf+1, lnh+1 ) );
		B1 *= sqrt(gsl_sf_fact(lnc));
		B1 *= sqrt(gsl_sf_fact(lnd));

		//boost::multiprecision::float128 B2 = 0;
		long double B2 = 0;
		int kmin = max(0, p-q-nb);
		int kmax = min(na, p-q);
		for (int k=kmin;k<=kmax;++k)
		{
			long double temp = 1;
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
			temp *= boost::math::tgamma_ratio( la+k +1, 2*la+2*k+1 +1);
			temp *= boost::math::tgamma_ratio( p-int((la-lb)/2)-k +1, 2*p-la+lb-2*k+1 +1);
			B2 += temp; 
		}
		double retB = (B1 * B2);//.convert_to<double>();
		return retB; //B1 * B2;
	}

	double
	RadialIntegral_RpowK(int na, int la, int nb, int lb, int k)
	{
		double I = 0;
		int pmin = (la+lb)/2;
		int pmax = pmin + na + nb;
		for (int p=pmin;p<=pmax;++p)
		{
			I += TalmiB(na,la,nb,lb,p) * TalmiI(p,k);
		}
		return I;
	} 

	int
	main()
	{
		int na = 0;
		int la = 0;
		int nb = 0;
		int lb = 0;
		int k = -1;
		double rad = RadialIntegral_RpowK( na,la, nb,lb, k);
		cout << "Rad =" << rad << endl;
		return 0;
	}
}
