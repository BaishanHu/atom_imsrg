#ifndef imsrg_util_hh
#define imsrg_util_hh 1

#include "ModelSpace.hh"
#include "Operator.hh"
#include "HartreeFock.hh"
#include "IMSRGSolver.hh"
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_laguerre.h>
#include <gsl/gsl_sf_gamma.h>
#include <vector>

#define HBARC 197.3269718 // hc in MeV * fm
#define M_NUCLEON 938.9185 // average nucleon mass in MeV
#define M_ELECTRON 0.510998910 // electron mass in MeV; scale to eV?
#define BOHR_RADIUS 0.0529 // Bohr Radius in nm

namespace imsrg_util
{
 Operator OperatorFromString(ModelSpace& modelspace, string str);
 map<index_t,double> GetSecondOrderOccupations(Operator& H, int emax);

 Operator NumberOp(ModelSpace& modelspace, int n, int l, int j2, int tz2);
 Operator NumberOpAlln(ModelSpace& modelspace, int l, int j2, int tz2);
 Operator PSquaredOp(ModelSpace& modelspace);
 Operator RSquaredOp(ModelSpace& modelspace);
 Operator E0Op(ModelSpace& modelspace);
 Operator ElectricMultipoleOp(ModelSpace& modelspace, int L);
 Operator MagneticMultipoleOp(ModelSpace& modelspace, int L);
 Operator MagneticMultipoleOp_pn(ModelSpace& modelspace, int L, string pn);
 Operator Trel_Op(ModelSpace& modelspace);
 Operator TCM_Op(ModelSpace& modelspace);
 Operator HCM_Op(ModelSpace& modelspace);
 Operator HarmonicOneBody(ModelSpace& modelspace);
 Operator InverseR_Op(ModelSpace& modelspace);
 Operator KineticEnergy_Op(ModelSpace& modelspace);
 Operator Energy_Op(ModelSpace& modelspace);

 Operator R2CM_Op(ModelSpace& modelspace);
 Operator Rp2_corrected_Op(ModelSpace& modelspace, int A, int Z);
 Operator Rn2_corrected_Op(ModelSpace& modelspace, int A, int Z);
 Operator Rm2_corrected_Op(ModelSpace& modelspace, int A, int Z);
 Operator R2_p1_Op(ModelSpace& modelspace);
 Operator R2_1body_Op(ModelSpace& modelspace, string option);
 Operator R2_p2_Op(ModelSpace& modelspace);
 Operator R2_2body_Op(ModelSpace& modelspace, string option);
 Operator ProtonDensityAtR(ModelSpace& modelspace, double R);
 Operator NeutronDensityAtR(ModelSpace& modelspace, double R);
 Operator RpSpinOrbitCorrection(ModelSpace& modelspace);
 Operator FourierBesselCoeff(ModelSpace& modelspace, int nu, double R, vector<index_t> index_list);

 Operator NumericalE2b(ModelSpace& modelspace);
 Operator CorrE2b(ModelSpace& modelspace);
 Operator CorrE2b_Hydrogen(ModelSpace& modelspace);
 double Corr_Invr(ModelSpace& modelspace, Ket & bra, Ket & ket, int J, string systemBasis);
 double Corr_Invr_Hydrogen(ModelSpace& modelspace, Ket & bra, Ket & ket, int J);

 void GenerateRadialIntegrals(ModelSpace& modelspace, int ind);
 unsigned long long int getMoshkey( int N, int Lam, int n, int lam, int n1, int l1, int n2, int l2, int L );
 unsigned long long int getNineJkey(double j1, double j2, double J12, double j3, double j4, double J34, double J13, double J24, double J);
 void PrecalculateRad_fromList( vector<unsigned int>& rad_list, ModelSpace& modelspace);
 double getRadialIntegral(int n1, int l1, int n2, int l2, ModelSpace& modelspace);


 Operator Isospin2_Op(ModelSpace& modelspace);
 Operator AllowedFermi_Op(ModelSpace& modelspace);
 Operator AllowedGamowTeller_Op(ModelSpace& modelspace);
 Operator Sigma_Op(ModelSpace& modelspace);
 Operator Sigma_Op_pn(ModelSpace& modelspace, string pn);
 Operator RadialOverlap(ModelSpace& modelspace);
 Operator LdotS_Op(ModelSpace& modelspace);
 unsigned long CalcCacheIndex(int J, int n1, int l1, int n2, int l2, int n3, int l3, int n4, int l4);
 unsigned long CalcLesserIndex(int n1, int l1, int n2, int l2);
 double Yml (double theta, int l, int m);
 double hydrogenWF(double x, double theta, int n, int l, int m, int Z, int A);
 int f(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);

 Operator eeCoulomb(ModelSpace& modelspace);
 Operator ElectronTwoBody(ModelSpace& modelspace);
 double CalculateCMInvR( double n1, double l1, double s1, double j1,
			 double n2, double l2, double s2, double j2,
			 double n3, double l3, double s3, double j3,
			 double n4, double l4, double s4, double j4, ModelSpace& modelspace, double J);

 Operator Single_Ref_1B_Density_Matrix(ModelSpace& modelspace); // This doesn't work
 double Get_Charge_Density(Operator& DM, double r);  // This doesn't work

 double Calculate_p1p2(ModelSpace& modelspace, Ket & bra, Ket & ket, int J);
 void Calculate_p1p2_all(Operator& OpIn);
 double Calculate_r1r2(ModelSpace& modelspace, Ket & bra, Ket & ket, int J);
 double HO_density(int n, int l, double hw, double r);
 double HO_Radial_psi(int n, int l, double hw, double r);
 double RadialIntegral(int na, int la, int nb, int lb, int L, ModelSpace& ms);
 double RadialIntegral_RpowK(int na, int la, int nb, int lb, int k, ModelSpace& ms);
 double TalmiI(int p, double k);
 double TalmiB(int na, int la, int nb, int lb, int p, ModelSpace& ms);
 vector<double> GetOccupationsHF(HartreeFock& hf);
 vector<double> GetOccupations(HartreeFock& hf, IMSRGSolver& imsrgsolver);
 vector<double> GetDensity(vector<double>& occ, vector<double>& R, vector<int>& orbits, ModelSpace& modelspace);

 void Embed1BodyIn2Body(Operator& op1, int A);
 double GetEmbeddedTBME(Operator& op1, index_t i, index_t j, index_t k, index_t l, int Jbra,int Jket, int Lambda);

 void CommutatorTest(Operator& X, Operator& Y);
 void Reduce(Operator&);
 void UnReduce(Operator&);
 /*
 extern "C"
 {
    #include "cube_test.h"
    //int fcube( int n1,int l1,int j1, int n2,int l2,int j2, int n3,int l3,int j3, int n4,int l4,int j4, int Z);
 } */

 //double Stirling(double n){return sqrt( 2 * 3.14159265359 * n ) * pow(n,n) / exp(n);}; // Stirling's Approximation
 vector<double> ser(double j); // returns a vector which is {1,2,...,j}
 bool isOnes(vector<double> a); // checks to see if a vector contains only ones.
 double simplefact(vector<double> n, vector<double> d, bool isSquare=false); // Reduces factorial division;
 //boost::multiprecision::cpp_bin_float_100 simplefact(vector<double> n, vector<double> d, bool isSquare=false); // Reduces factorial division;


// Templated functions need to be defined in the header file (or else explicitly declared in the .cc file).
 template <typename T>
 T VectorUnion(T& v1)
 {
   return v1;
 }
 
 template <typename T, typename... Args>
 T VectorUnion(T& v1, T& v2, Args... args)
 {
   T vec(v1.size()+v2.size());
   copy(v1.begin(),v1.end(),vec.begin());
   copy(v2.begin(),v2.end(),vec.begin()+v1.size());
   return VectorUnion(vec, args...);
 }

}




#endif
