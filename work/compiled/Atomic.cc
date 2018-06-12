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
  
  cout << "SystemBasis: " << systemBasis << endl;

  // test 2bme file
  ifstream test(inputtbme);
  if( not test.good() )
  {
    cout << "trouble reading " << inputtbme << " exiting. " << endl;
    //return 1;
  }
  test.close();
  // test 3bme file
  //if (input3bme != "none")
  //{
  //  test.open(input3bme);
  //  if( not test.good() )
  //  {
  //    cout << "trouble reading " << input3bme << " exiting. " << endl;
  //    return 1;
  //  }
  //  test.close();
  //}


  //if (Lmax == 0) {
  //  cout << "Lmax not recognized, setting to Lmax=2" << endl;
  //  Lmax = 2;
  //}
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

  //for (auto& orb : modelspace.Orbits) {
  //  cout << "orb.index=" << orb.index << " orb.n=" << orb.n << " orb.l=" << orb.l << " orb.j2=" << orb.j2 << " orb.occ=" << orb.occ << endl;
  //}

  cout << "Setting particle_rank." << endl;
  int particle_rank = input3bme=="none" ? 2 : 3;
  cout << "particle_rank=" << particle_rank << " Constructing Hbare, norbits=" << modelspace.GetNumberOrbits() << endl;
  //Operator Hbare = Operator(modelspace,0,0,0,particle_rank);
  Operator Hbare = Operator(modelspace);
  cout << "Constructed Hbare, setting hermitian." << endl;
  Hbare.SetHermitian();
  cout << "Set Hermitian, destroying any relic 0, 1, or 2 body elements." << endl;
  Hbare.Erase();
  cout << "Erased relic 0, 1, and 2 body elements about to ReadBareTBME_Darmstadt" << endl;
  if (use_brueckner_bch == "true" or use_brueckner_bch == "True")
  {
    Hbare.SetUseBruecknerBCH(true);
    cout << "Using Brueckner flavor of BCH" << endl;
  }

//  cout << "Reading interactions..." << endl;
//
//  #pragma omp parallel sections 
//  {
//    #pragma omp section
//    {
//    if (fmt2 == "me2j")
//      rw.ReadBareTBME_Darmstadt(inputtbme, Hbare,file2e1max,file2e2max,file2lmax);
//    else if (fmt2 == "navratil" or fmt2 == "Navratil")
//      rw.ReadBareTBME_Navratil(inputtbme, Hbare);
//    else if (fmt2 == "oslo" )
//      rw.ReadTBME_Oslo(inputtbme, Hbare);
//    else if (fmt2 == "oakridge" )
//      rw.ReadTBME_OakRidge(inputtbme, Hbare);
//     cout << "done reading 2N" << endl;
//    }
//  
//    #pragma omp section
//    if (Hbare.particle_rank >=3)
//    {
//      rw.Read_Darmstadt_3body(input3bme, Hbare, file3e1max,file3e2max,file3e3max);
//      cout << "done reading 3N" << endl;
//    }  
//  }
  //cout << "Modelspace has this many tbc: " << modelspace.nTwoBodyChannels << endl;
  //rw.ReadBareTBME_Darmstadt( inputtbme, Hbare, file2e1max, file2e2max, file2lmax);
  
  if (systemBasis == "harmonic") {
    cout << "About to precalculate factorials for m=2*(2*emax + lmax)=" << 2*(2*eMax + 1*Lmax) << endl;
    modelspace.GenerateFactorialList( 2*(2*eMax + 1*Lmax));
    cout << "About to precalculate radial integrals." << endl;
    GenerateRadialIntegrals(modelspace,1*eMax*101011);
    cout << "FactorialList calculated." << endl;
    //rw.ReadOperatorFromJSON( inputtbme, Hbare, eMax, 2*eMax, Lmax, 1 );
    cout << "Read in interaction, moving to precalculating moshinsky." << endl;
    cout << "Adding T and V to Hbare; emax=" << modelspace.GetEmax() << endl;
    //Hbare += Trel_Op(modelspace) + InverseR_Op(modelspace);
    //Operator KE = KineticEnergy_Op(modelspace);
    //Operator& K = KE;
    //Operator IR = InverseR_Op(modelspace);
    //Operator& I = IR;
    //Hbare += KE + IR;
    cout << "TargetZ=" << modelspace.GetTargetZ() << endl;
    cout << "Adding KE to Hbare." << endl;
    Hbare += KineticEnergy_Op(modelspace);
    cout << "Added KE." << endl;
    cout << "Adding InvR to Hbare." << endl;
    Hbare += InverseR_Op(modelspace) * modelspace.GetTargetZ();
    cout << "Added InvR; adding two body." << endl;
    //Hbare += CorrE2b(modelspace);
    //cout << "Added Twobody, moving on." << endl;
  } else {
    //cout << "About to precalculate factorials for m=2*(2*emax + lmax)=" << 2*(2*eMax + 1*Lmax) << endl;
    //modelspace.GenerateFactorialList( 2*(2*eMax + 1*Lmax)+40 );
    //cout << "About to precalculate radial integrals." << endl;
    //GenerateRadialIntegrals(modelspace,2*eMax*101011);
    //cout << "FactorialList calculated." << endl;
    //modelspace.PreCalculateMoshinsky( systemBasis );
    //cout << "Precalculated Mosh, moving on." << endl;
    //rw.ReadOperatorFromJSON( inputtbme, Hbare, eMax, 2*eMax, Lmax, 1 );
    //modelspace.GenerateOsToHydroCoeff(eMax);
    //Hbare += NumericalE2b(modelspace);
    //Hbare += CorrE2b_Hydrogen(modelspace);
    cout << "Adding Hydrogen Energies." << endl;
    Hbare += Energy_Op(modelspace);
    cout << "Onebody:" << endl << Hbare.OneBody << endl;
    cout << "Adding two-body correction." << endl;
    Hbare += ElectronTwoBody(modelspace);
  }
/*
  cout << "OneBody=" << endl << Hbare.OneBody << endl;
  cout << "TwoBody=" << endl;
  for (int ch = 0; ch < Hbare.nChannels; ch++) {
    cout << "----- Channel " << ch << " with J=" << modelspace.GetTwoBodyChannel(ch).J << "-----" << endl;
    Hbare.PrintTwoBody(ch);
    cout << endl;
  }
*/

  //cout << "Adding ElectronTwoBody to Hbare." << endl;
  //Hbare += CorrE2b(modelspace);
  //Hbare += ElectronTwoBody(modelspace);
  //cout << "Added ElectronTwoBody to Hbare." << endl;
  //cout << "OneBody=" << Hbare.OneBody << endl;
  //for (auto& i : K.OneBodyChannels)
  //{
     
  //}

  if (abs(BetaCM) > 1e-3)
  {
    Hbare += BetaCM * HCM_Op(modelspace);
  }
  //Hbare.OneBody.print();
  cout << "About to create hf(Hbare)" << endl;
  HartreeFock hf(Hbare);
  cout << "done Converting Hbare to HF basis" << endl;
  hf.Solve();
  cout << "EHF = " << hf.EHF << endl;

//  Hbare -= BetaCM * 1.5*hw;

  if (method != "HF")
  {
    cout << "Perturbative estimates of gs energy:" << endl;
    double EMP2 = Hbare.GetMP2_Energy();
    cout << "EMP2 = " << EMP2 << endl; 
    double EMP3 = Hbare.GetMP3_Energy();
    cout << "EMP3 = " << EMP3 << endl; 
    cout << "To 3rd order, E = " << Hbare.ZeroBody+EMP2+EMP3 << endl;
  }
  cout << "About to calculate all of the operators." << endl;
  // Calculate all the desired operators
  for (auto& opname : opnames)
  {
           if (opname == "R2_p1")        ops.emplace_back( R2_1body_Op(modelspace,"proton") );
      else if (opname == "R2_p2")        ops.emplace_back( R2_2body_Op(modelspace,"proton") );
      else if (opname == "R2_n1")        ops.emplace_back( R2_1body_Op(modelspace,"neutron") );
      else if (opname == "R2_n2")        ops.emplace_back( R2_2body_Op(modelspace,"neutron") );
      else if (opname == "Rp2")          ops.emplace_back( Rp2_corrected_Op(modelspace,modelspace.GetTargetMass(),modelspace.GetTargetZ()) );
      else if (opname == "Rn2")          ops.emplace_back( Rn2_corrected_Op(modelspace,modelspace.GetTargetMass(),modelspace.GetTargetZ()) );
      else if (opname == "Rm2")          ops.emplace_back( Rm2_corrected_Op(modelspace,modelspace.GetTargetMass(),modelspace.GetTargetZ()) );
      else if (opname == "TCM_Op")     	 ops.emplace_back( TCM_Op(modelspace) );
      else if (opname == "KineticEnergy")ops.emplace_back( KineticEnergy_Op(modelspace) );
      else if (opname == "InverseR")     ops.emplace_back( InverseR_Op(modelspace) );
      else if (opname == "ElectronTwoBody") ops.emplace_back( ElectronTwoBody(modelspace) );
      else if (opname == "CorrE2b")	 ops.emplace_back( CorrE2b(modelspace) );
      else if (opname == "CorrE2b_Hydrogen")	 ops.emplace_back( CorrE2b_Hydrogen(modelspace) );
      //else if (opname == "NumericalE2b") ops.emplace_back( NumericalE2b(modelspace) );
      else if (opname == "Energy_Op")    ops.emplace_back( Energy_Op(modelspace) );
      else if (opname == "E2")           ops.emplace_back( ElectricMultipoleOp(modelspace,2) );
      else if (opname == "M1")           ops.emplace_back( MagneticMultipoleOp(modelspace,1) );
      else if (opname == "Fermi")        ops.emplace_back( AllowedFermi_Op(modelspace) );
      else if (opname == "GamowTeller")  ops.emplace_back( AllowedGamowTeller_Op(modelspace) );
      else if (opname == "R2CM")         ops.emplace_back( R2CM_Op(modelspace) );
      else if (opname == "HCM")          ops.emplace_back( HCM_Op(modelspace) );
      else if (opname == "Rso")          ops.emplace_back( RpSpinOrbitCorrection(modelspace) );
      else if (opname.substr(0,4) == "HCM_") // GetHCM with a different frequency, ie HCM_24 for hw=24
      {
         double hw_HCM;
//         double hw_save = modelspace.GetHbarOmega();
         istringstream(opname.substr(4,opname.size())) >> hw_HCM;
//         modelspace.SetHbarOmega(hw_HCM);
//         ops.emplace_back( HCM_Op(modelspace) );
         int A = modelspace.GetTargetMass();
         Operator hcm = TCM_Op(modelspace) + 0.5*A*M_NUCLEON*hw*hw/HBARC/HBARC*R2CM_Op(modelspace); 
         ops.emplace_back( hcm );
//         modelspace.SetHbarOmega(hw_save);
      }
      else if (opname.substr(0,4) == "Rp2Z")
      {
        int Z_rp;
        istringstream(opname.substr(4,opname.size())) >> Z_rp;
        ops.emplace_back( Rp2_corrected_Op(modelspace,modelspace.GetTargetMass(),Z_rp) );
      }
      else if (opname.substr(0,4) == "Rn2Z")
      {
        int Z_rp;
        istringstream(opname.substr(4,opname.size())) >> Z_rp;
        ops.emplace_back( Rn2_corrected_Op(modelspace,modelspace.GetTargetMass(),Z_rp) );
      }
      else if (opname.substr(0,4) == "rhop")
      {
        double rr;
        istringstream(opname.substr(4,opname.size())) >> rr;
        ops.emplace_back( ProtonDensityAtR(modelspace,rr));
      }
      else if (opname.substr(0,4) == "rhon")
      {
        double rr;
        istringstream(opname.substr(4,opname.size())) >> rr;
        ops.emplace_back( NeutronDensityAtR(modelspace,rr));
      }
      else if (opname.substr(0,6) == "OneOcc")
      {
         map<char,int> lvals = {{'s',0},{'p',1},{'d',2},{'f',3},{'g',4},{'h',5}};
         char pn,lspec;
         int n,l,j,t;
         istringstream(opname.substr(6,1)) >> pn;
         istringstream(opname.substr(7,1)) >> n;
         istringstream(opname.substr(8,1)) >> lspec;
         istringstream(opname.substr(9,opname.size())) >> j;
         l = lvals[lspec];
         t = pn == 'p' ? -1 : 1;
         ops.emplace_back( NumberOp(modelspace,n,l,j,t) );
      }
      else if (opname.substr(0,6) == "AllOcc")
      {
         map<char,int> lvals = {{'s',0},{'p',1},{'d',2},{'f',3},{'g',4},{'h',5}};
         char pn,lspec;
         int l,j,t;
         istringstream(opname.substr(6,1)) >> pn;
         istringstream(opname.substr(7,1)) >> lspec;
         istringstream(opname.substr(8,opname.size())) >> j;
         l = lvals[lspec];
         t = pn == 'p' ? -1 : 1;
         ops.emplace_back( NumberOpAlln(modelspace,l,j,t) );
      }
      else if (opname.substr(0,9) == "protonFBC")
      {
         int nu;
         istringstream(opname.substr(9,opname.size())) >> nu;
         ops.emplace_back( FourierBesselCoeff( modelspace, nu, 8.0, modelspace.proton_orbits) );
      }
      else if (opname.substr(0,10) == "neutronFBC")
      {
         int nu;
         istringstream(opname.substr(10,opname.size())) >> nu;
         ops.emplace_back( FourierBesselCoeff( modelspace, nu, 8.0, modelspace.neutron_orbits) );
      }
      else //need to remove from the list
      {
         cout << "Unknown operator: " << opname << endl;
      }
  }

  for (int i=0;i<ops.size();++i)
  {
    Operator& op = ops[i];
    cout << opnames[i] << " has OneBody matrix (before diag):" << endl;
    cout << ops[i].OneBody << endl;
    //cout << opnames[i] << " has TwoBody matrix (before diag):" << endl;
    //cout << ops[i].TwoBody << endl;
  }
  cout << "In Atomic, about to normal order Hbare." << endl;
  if (basis == "HF" and method !="HF")
    Hbare = hf.GetNormalOrderedH();
  else if (basis == "oscillator")
    Hbare = Hbare.DoNormalOrdering();
  cout << "Hbare should be normal ordered." << endl;
  

  for (auto& op : ops)
  {
     if (basis == "HF") op = hf.TransformToHFBasis(op);
     op = op.DoNormalOrdering();
     if (method == "MP3")
     {
       double dop = op.MP1_Eval( Hbare );
       cout << "Operator 1st order correction  " << dop << "  ->  " << op.ZeroBody + dop << endl;
     }
  }
  auto itR2p = find(opnames.begin(),opnames.end(),"Rp2");
  if (itR2p != opnames.end())
  {
    Operator& Rp2 = ops[itR2p-opnames.begin()];
    int Z = modelspace.GetTargetZ();
    int A = modelspace.GetTargetMass();
    cout << " HF point proton radius = " << sqrt( Rp2.ZeroBody ) << endl; 
    cout << " HF charge radius = " << sqrt( Rp2.ZeroBody + r2p + r2n*(A-Z)/Z + DF) << endl; 
  }
  for (int i=0;i<ops.size();++i)
  {
    Operator& op = ops[i];
    cout << opnames[i] << " = " << ops[i].ZeroBody << endl;
    cout << opnames[i] << " has full matrix (after transforming):" << endl;
    cout << "One Body = " << endl;
    cout << ops[i].OneBody << endl;
    //cout << "Two Body = " << endl;
  }
  
  if ( method == "HF" or method == "MP3")
  {
    hf.PrintSPE();
    Hbare.PrintTimes();
    return 0;
  }
  cout << "About to set IMSRGSolver for Hbare." << endl;
  IMSRGSolver imsrgsolver(Hbare);
  imsrgsolver.SetReadWrite(rw);
  
  if (method == "NSmagnus") // "No split" magnus
  {
    omega_norm_max=500;
    method = "magnus";
  }
  imsrgsolver.SetMethod(method);
  imsrgsolver.SetHin(Hbare);
  imsrgsolver.SetSmax(smax);
  imsrgsolver.SetFlowFile(flowfile);
  imsrgsolver.SetDs(ds_0);
  imsrgsolver.SetDenominatorDelta(denominator_delta);
  imsrgsolver.SetdOmega(domega);
  imsrgsolver.SetOmegaNormMax(omega_norm_max);
  imsrgsolver.SetODETolerance(ode_tolerance);
  if (denominator_delta_orbit != "none")
    imsrgsolver.SetDenominatorDeltaOrbit(denominator_delta_orbit);

  if (nsteps > 1) // two-step decoupling, do core first
  {
    imsrgsolver.SetGenerator(core_generator);
    imsrgsolver.Solve();
    if (method == "magnus") smax *= 2;
  }
  cout << "About to Set valence_generator." << endl;
  imsrgsolver.SetGenerator(valence_generator);
  cout << "About to set smax." << endl;
  imsrgsolver.SetSmax(smax);
  cout << "About to Solve imsrg." << endl;
  imsrgsolver.Solve();


  // Transform all the operators
  if (method == "magnus")
  {
    if (ops.size()>0) cout << "transforming operators" << endl;
    for (size_t i=0;i<ops.size();++i)
    {
      cout << opnames[i] << " " << flush;
      ops[i] = imsrgsolver.Transform(ops[i]);
      cout << " (" << ops[i].ZeroBody << " ) " << endl; 
    }
    cout << endl;
    // increase smax in case we need to do additional steps
    smax *= 1.5;
    imsrgsolver.SetSmax(smax);
  }


  // If we're doing targeted normal ordering 
  // we now re-normal order wrt to the core
  // and do any remaining flow.
//  if (reference != "default"  and reference != valence_space)
  ModelSpace ms2(modelspace);
  bool renormal_order = false;
  if (modelspace.valence.size() > 0 )
  {
    renormal_order = modelspace.holes.size() != modelspace.core.size();
    if (not renormal_order)
    {
      for (auto c : modelspace.core)
      {
         if ( (find( modelspace.holes.begin(), modelspace.holes.end(), c) == modelspace.holes.end()) or (abs(1-modelspace.holes[c])>1e-6))
         {
           renormal_order = true;
           break;
         }
      }
    }
  }
//  if ( modelspace.core != modelspace.holes )
  if ( renormal_order )
  {

    Hbare = imsrgsolver.GetH_s();

    int nOmega = imsrgsolver.GetOmegaSize() + imsrgsolver.GetNOmegaWritten();
    cout << "Undoing NO wrt A=" << modelspace.GetAref() << " Z=" << modelspace.GetZref() << endl;
    Hbare = Hbare.UndoNormalOrdering();

    ms2.SetReference(ms2.core); // chage the reference determinant
    Hbare.SetModelSpace(ms2);

    cout << "Doing NO wrt A=" << ms2.GetAref() << " Z=" << ms2.GetZref() << "  norbits = " << ms2.GetNumberOrbits() << endl;
    Hbare = Hbare.DoNormalOrdering();

    imsrgsolver.SetHin(Hbare);
    imsrgsolver.SetEtaCriterion(1e-4);
    imsrgsolver.Solve();
    // Change operators to the new basis, then apply the rest of the transformation
    cout << "Final transformation on the operators..." << endl;
    for (auto& op : ops)
    {
      double ZeroBody_before = op.ZeroBody;
      op = op.UndoNormalOrdering();
      double ZeroBody_undo = op.ZeroBody;
      op.SetModelSpace(ms2);
      op = op.DoNormalOrdering();
      double ZeroBody_mid = op.ZeroBody;
      // transform using the remaining omegas
      op = imsrgsolver.Transform_Partial(op,nOmega);
      cout << ZeroBody_before << "   =>   " << ZeroBody_undo << "   =>   " << ZeroBody_mid<< "   =>   " << op.ZeroBody << endl;
    }
  }


  // Write the output

  // If we're doing a shell model interaction, write the
  // interaction files to disk.
  if (modelspace.valence.size() > 0)
  {
    if (valence_file_format == "antoine")
    {
      rw.WriteAntoine_int(imsrgsolver.GetH_s(),intfile+".ant");
      rw.WriteAntoine_input(imsrgsolver.GetH_s(),intfile+".inp",modelspace.GetAref(),modelspace.GetZref());
    }
//    else
//    {
      rw.WriteNuShellX_int(imsrgsolver.GetH_s(),intfile+".int");
      rw.WriteNuShellX_sps(imsrgsolver.GetH_s(),intfile+".sp");
//    }

    if (method == "magnus")
    {
       for (int i=0;i<ops.size();++i)
       {
//          ops[i] = imsrgsolver.Transform(ops[i]);
          if ((ops[i].GetJRank()+ops[i].GetTRank()+ops[i].GetParity())<1)
          {
            rw.WriteNuShellX_op(ops[i],intfile+opnames[i]+".int");
          }
          else
          {
            rw.WriteTensorOneBody(intfile+opnames[i]+"_1b.op",ops[i],opnames[i]);
            rw.WriteTensorTwoBody(intfile+opnames[i]+"_2b.op",ops[i],opnames[i]);
          }
       }
    }
  }
  else // single ref. just print the zero body pieces out. (maybe check if its magnus?)
  {
    cout << "Core Energy = " << setprecision(6) << imsrgsolver.GetH_s().ZeroBody << endl;
    for (int i=0;i<ops.size();++i)
    {
      Operator& op = ops[i];
      cout << opnames[i] << " = " << ops[i].ZeroBody << endl;
      if ( opnames[i] == "Rp2" )
      {
         int Z = modelspace.GetTargetZ();
         int A = modelspace.GetTargetMass();
         cout << " IMSRG point proton radius = " << sqrt( op.ZeroBody ) << endl; 
         cout << " IMSRG charge radius = " << sqrt( op.ZeroBody + r2p + r2n*(A-Z)/Z + DF) << endl; 
      }
    }
  }


 
 
  Hbare.PrintTimes();
  cout << "That's all, folks!" << endl;
  return 0;

}
