#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <omp.h>
#include "IMSRG.hh"
#include "Parameters.hh"

using namespace imsrg_util;

int main(int argc, char** argv)
{
  // Default parameters, and everything passed by command line args.
  Parameters PAR(argc,argv);

  string inputtbme = PAR.s("2bme");
  string input3bme = PAR.s("3bme");
  string reference = PAR.s("reference");
  string valence_space = PAR.s("valence_space");
  string basis = PAR.s("basis");
  string systemBasis = PAR.s("systemBasis");
  string method = PAR.s("method");
  string flowfile = PAR.s("flowfile");
  string intfile = PAR.s("intfile");
  string core_generator = PAR.s("core_generator");
  string valence_generator = PAR.s("valence_generator");
  string fmt2 = PAR.s("fmt2");
  string denominator_delta_orbit = PAR.s("denominator_delta_orbit");
  string LECs = PAR.s("LECs");
  string scratch = PAR.s("scratch");
  string use_brueckner_bch = PAR.s("use_brueckner_bch");
  string valence_file_format = PAR.s("valence_file_format");
  string systemtype = PAR.s("systemtype");

  int eMax = PAR.i("emax");
  int Lmax = PAR.i("Lmax");
  int E3max = PAR.i("e3max");
  int lmax3 = PAR.i("lmax3");
  int targetMass = PAR.i("A");
  int nsteps = PAR.i("nsteps");
  int file2e1max = PAR.i("file2e1max");
  int file2e2max = PAR.i("file2e2max");
  int file2lmax = PAR.i("file2lmax");
  int file3e1max = PAR.i("file3e1max");
  int file3e2max = PAR.i("file3e2max");
  int file3e3max = PAR.i("file3e3max");

  double hw = PAR.d("hw");
  double smax = PAR.d("smax");
  double ode_tolerance = PAR.d("ode_tolerance");
  double ds_max = PAR.d("ds_max");
  double ds_0 = PAR.d("ds_0");
  double domega = PAR.d("domega");
  double omega_norm_max = PAR.d("omega_norm_max"); 
  double denominator_delta = PAR.d("denominator_delta");
  double BetaCM = PAR.d("BetaCM");

  vector<string> opnames = PAR.v("Operators");

  vector<Operator> ops;

  ReadWrite rw;
  rw.SetLECs_preset(LECs);
  rw.SetScratchDir(scratch);
  //string SystemType = systemtype==string::empty "atomic" : systemtype;

  string SystemType = "atomic";
  if (Lmax > eMax) Lmax = eMax;
  cout << "About to construct modelspace eMax="<< eMax << " Lmax=" << Lmax << " SystemType=" << SystemType << " Valence_Space=" << valence_space << " reference=" << reference << " systemBasis=" << systemBasis << endl;
  if (eMax == 0) eMax = 4;
  if (reference == "default") reference = "He2";
  //if (valence_space == "016")
  ModelSpace modelspace = ModelSpace(eMax, reference, valence_space, Lmax, SystemType, systemBasis);
  cout << "Default modelspace constructed, setting SystemType." << endl;
  modelspace.SetSystemType(SystemType);
  modelspace.SetSystemBasis(systemBasis);
  //cout << "About to precalculate factorials for m=2*(2*emax + lmax)=" << 2*(2*eMax + 1*Lmax) << endl;
  //modelspace.GenerateFactorialList( 2*(2*eMax + 1*Lmax) );
  //cout << "FactorialList calculated." << endl;
  //cout << "About to generate radial integrals." << endl;
  //GenerateRadialIntegrals(modelspace,2*eMax*101011);
  //cout << "Generated Radial integrals." << endl;
  cout << "About to generate factorialList." << endl;
  modelspace.GenerateFactorialList( min(4*(2*eMax + 1*Lmax),170) );
  cout << "About to precalculate radial integrals." << endl;
  GenerateRadialIntegrals(modelspace, eMax*1010000 + (eMax+Lmax)*202); // Should handle most/all integrals

  cout << "modelspace initialized." << endl;
  cout << "Emax is "<< modelspace.GetEmax() << endl;

  if (nsteps < 0)
    nsteps = modelspace.valence.size()>0 ? 2 : 1;


  modelspace.SetHbarOmega(hw);
  cout << "Setting target mass." << endl;
  //int targetMass = 0;
  if (targetMass>0)
     modelspace.SetTargetMass(targetMass);
  cout << "Setting E3Max." << endl;
  modelspace.SetE3max(E3max);
  cout << "Setting lmax3." << endl;
  if (lmax3>0)
     modelspace.SetLmax3(lmax3);

  cout << "Setting particle_rank." << endl;
  int particle_rank = input3bme=="none" ? 2 : 3;
  cout << "particle_rank=" << particle_rank << " Constructing Hbare, norbits=" << modelspace.GetNumberOrbits() << endl;
  //Operator Hbare = Operator(modelspace,0,0,0,particle_rank);
  Operator Hbare = Operator(modelspace);
  cout << "Constructed Hbare, setting hermitian." << endl;
  Hbare.SetHermitian();
  cout << "Set Hermitian, destroying any relic 0, 1, or 2 body elements." << endl;
  Hbare.Erase();

  if (use_brueckner_bch == "true" or use_brueckner_bch == "True")
  {
    Hbare.SetUseBruecknerBCH(true);
    cout << "Using Brueckner flavor of BCH" << endl;
  }

  //cout << "Added InvR." << endl;
  //modelspace.PreCalculateMoshinsky();
  //cout << "Precalculated Mosh, moving on." << endl;
  cout << "Adding ElectronTwoBody to Hbare." << endl;
  if ( systemBasis == "hydrogen" ) {
    //modelspace.GenerateOsToHydroCoeff(eMax);
    Hbare += ElectronTwoBody(modelspace);
    //Hbare += CorrE2b_Hydrogen(modelspace);
  } else Hbare += CorrE2b(modelspace);
  //Hbare += ElectronTwoBody(modelspace);
  cout << "Added ElectronTwoBody to Hbare." << endl;
  std::stringstream fn;
  //Operator New = Operator(modelspace);
  if (systemBasis == "harmonic") {
	fn << "/home/dlivermore/ragnar_imsrg/work/scripts/atomicME_" << reference << "_basis_" << systemBasis << "Aug30_emax" << eMax << "_hw" << hw << ".me2j";
	cout << "Writing Hbare to file with filename=" << fn.str() << endl;
	rw.Write_me2j( fn.str(),  Hbare, eMax, 3*eMax, -1 );
	/*
	for (int ch = 0; ch < Hbare.nChannels; ch++) {
        	cout << endl;
        	cout << "----- Channel " << ch << " -----" << endl;
        	Hbare.PrintTwoBody(ch);
        	cout << endl;
  	} */

	cout << "Written, making new operator." << endl;
	//rw.ReadBareTBME_Darmstadt( fn.str(), New, eMax, 3*eMax, -1 );
  } else {
  	fn << "/home/dlivermore/ragnar_imsrg/work/scripts/atomicME_" << reference << "_basis_" << systemBasis << "Aug30_emax" << eMax << "_lmax" << Lmax << "_hw" << hw << ".json";
	cout << "Writing Hbare to file with filename=" << fn.str() << endl;
	rw.WriteOperatorToJSON( fn.str(), Hbare, eMax, 3*eMax, Lmax, 0.1 );
	cout << "Written, making new operator." << endl;
	cout << "Reading back in operator." << endl;
	//rw.ReadOperatorFromJSON( fn.str(), New, eMax, 3*eMax, Lmax, 0.1 );
  }

  //Operator Diff = Hbare - New;
  //cout << "Norm of New=" << New.Norm() << endl;
  cout << "Norm of Hbare=" << Hbare.Norm() << endl;
  //cout << "Norm of diff=" << Diff.Norm() << endl;

  cout << endl;
  //cout << "OneBodyNorm of New=" << New.OneBodyNorm() << endl;
  cout << "OneBodyNorm of Hbare=" << Hbare.OneBodyNorm() << endl;
  //cout << "OneBodyNorm of diff=" << Diff.OneBodyNorm() << endl;

  cout << endl;
  //cout << "TwoBodyNorm of New=" << New.TwoBodyNorm() << endl;
  cout << "TwoBodyNorm of Hbare=" << Hbare.TwoBodyNorm() << endl;
  //cout << "TwoBodyNorm of diff=" << Diff.TwoBodyNorm() << endl;
 /*
  cout << "New:" << endl;
  cout << "OneBody=" << endl << New.OneBody << endl;
  cout << "TwoBody=" << endl;
  for (int ch = 0; ch < New.nChannels; ch++) {
    cout << "----- Channel " << ch << " with J=" << modelspace.GetTwoBodyChannel(ch).J << "-----" << endl;
    New.PrintTwoBody(ch);
    cout << endl;
  }
  cout << endl << endl;

  cout << "Hbare:" << endl;
  cout << "OneBody=" << endl << Hbare.OneBody << endl;
  cout << "TwoBody=" << endl;
  for (int ch = 0; ch < Hbare.nChannels; ch++) {
    cout << "----- Channel " << ch << " with J=" << modelspace.GetTwoBodyChannel(ch).J << "-----" << endl;
    Hbare.PrintTwoBody(ch);
    cout << endl;
  } 

  cout << "Diff:" << endl;
  cout << "OneBody=" << endl << New.OneBody << endl;
  cout << "TwoBody=" << endl;
  for (int ch = 0; ch < Diff.nChannels; ch++) {
    cout << "----- Channel " << ch << " with J=" << modelspace.GetTwoBodyChannel(ch).J << "-----" << endl;
    Diff.PrintTwoBody(ch);
    cout << endl;
  }
  cout << endl << endl;
 */
 /* for (int ch = 0; ch < Hbare.nChannels; ch++) {
	cout << endl;
	cout << "----- Channel " << ch << " -----" << endl;
	cout << endl;
	cout << "Hbare:" << endl;
	Hbare.PrintTwoBody(ch);
	cout << endl << "New:" << endl;
	New.PrintTwoBody(ch);
	cout << endl << "Diff:" << endl;
	Diff.PrintTwoBody(ch);
  } */
  Hbare.PrintTimes();
  cout << endl << "That's all, folks!" << endl;
  return 0;

}
