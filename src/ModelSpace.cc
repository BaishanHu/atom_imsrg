#include "ModelSpace.hh"
#include "AngMom.hh"
#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <gsl/gsl_math.h>

#include "imsrg_util.hh"
#include "AngMom.hh"
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/implicit_cast.hpp>
#include <gsl/gsl_integration.h>
#include <list>

#include "ModelSpace.hh"
#include "Operator.hh"
#include "HartreeFock.hh"
#include "IMSRGSolver.hh"
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_laguerre.h>
#include <gsl/gsl_sf_gamma.h>
#include <vector>
#include <cmath>

using namespace std;

Orbit::~Orbit()
{
//  cout << "In Orbit destructor" << endl;
}

Orbit::Orbit()
: n(-1), l(-1), j2(-1), tz2(-1),occ(-1),cvq(-1),index(-1)
//: n(-1), l(-1), j2(-1), tz2(-1),ph(-1),io(-1),index(-1)
{}

Orbit::Orbit(int n, int l, int j2, int tz2, double occ, int cvq, int index)
: n(n), l(l), j2(j2), tz2(tz2),occ(occ),cvq(cvq),index(index)
//: n(n), l(l), j2(j2), tz2(tz2),ph(ph),io(io),index(index)
{}

Orbit::Orbit(const Orbit& orb)
: n(orb.n), l(orb.l), j2(orb.j2), tz2(orb.tz2),occ(orb.occ),cvq(orb.cvq),index(orb.index)
//: n(orb.n), l(orb.l), j2(orb.j2), tz2(orb.tz2),ph(orb.ph),io(orb.io),index(orb.index)
{}

//************************************************************************
//************************************************************************
//************************************************************************
Ket::~Ket()
{
//  cout << "In Ket destructor" << endl;
}

Ket::Ket()
{}

Ket::Ket(Orbit& op_in, Orbit& oq_in)
: op(&op_in), oq(&oq_in), p(op_in.index), q(oq_in.index)
{
   phase_prefactor = ((op->j2+oq->j2)/2 + 1) % 2==0 ? 1 : -1;
   dpq = p==q ? 1 : 0;
}

int Ket::Phase(int J)
{
   return phase_prefactor * (J%2==0 ? 1 : -1);
}

//************************************************************************
//************************************************************************
//************************************************************************

TwoBodyChannel::~TwoBodyChannel()
{
//  cout << "In TwoBodyChannel destructor" << endl;
}

TwoBodyChannel::TwoBodyChannel()
{}

TwoBodyChannel::TwoBodyChannel(int j, int p, int t, ModelSpace *ms)
{
   Initialize(ms->GetTwoBodyChannelIndex(j,p,t), ms);
}

TwoBodyChannel::TwoBodyChannel(int N, ModelSpace *ms)
{
   Initialize(N,ms);
}

void TwoBodyChannel::Initialize(int N, ModelSpace *ms)
{
   int tbjmax = ms->TwoBodyJmax;
   J = N%(tbjmax+1);
   parity = (N/(tbjmax+1))%2;
   Tz = (N/(2*(tbjmax+1))-1);
   modelspace = ms;
   NumberKets = 0;
   int nk = modelspace->GetNumberKets();
   KetMap.resize(nk,-1); // set all values to -1
   for (int i=0;i<nk;i++)
   {
      //cout << "Geting Ketmap init'd; i=" << i << endl;
      Ket &ket = modelspace->GetKet(i);
      if ( CheckChannel_ket(ket) )
      {
	 //cout << "Got into check_channel." << endl;
         KetMap[i] = NumberKets;
         KetList.push_back(i);
         NumberKets++;
      }
   }
   KetIndex_pp = GetKetIndexFromList(modelspace->KetIndex_pp);
   KetIndex_hh = GetKetIndexFromList(modelspace->KetIndex_hh);
   KetIndex_ph = GetKetIndexFromList(modelspace->KetIndex_ph);
   KetIndex_cc = GetKetIndexFromList(modelspace->KetIndex_cc);
   KetIndex_vc = GetKetIndexFromList(modelspace->KetIndex_vc);
   KetIndex_qc = GetKetIndexFromList(modelspace->KetIndex_qc);
   KetIndex_vv = GetKetIndexFromList(modelspace->KetIndex_vv);
   KetIndex_qv = GetKetIndexFromList(modelspace->KetIndex_qv);
   KetIndex_qq = GetKetIndexFromList(modelspace->KetIndex_qq);
   vector<double> occvec;
   vector<double> unoccvec;
   for (index_t i=0;i<modelspace->KetIndex_hh.size();++i)
   {
      if (CheckChannel_ket(modelspace->GetKet(modelspace->KetIndex_hh[i])))
      {
        occvec.push_back( modelspace->Ket_occ_hh[i]);
        unoccvec.push_back( modelspace->Ket_unocc_hh[i]);
      }
   }
   Ket_occ_hh = arma::vec(occvec);
   Ket_unocc_hh = arma::vec(unoccvec);
   occvec.clear();
   unoccvec.clear();
   for (index_t i=0;i<modelspace->KetIndex_ph.size();++i)
   {
      if (CheckChannel_ket(modelspace->GetKet(modelspace->KetIndex_ph[i])))
      {
        occvec.push_back( modelspace->Ket_occ_ph[i]);
        unoccvec.push_back( modelspace->Ket_unocc_ph[i]);
      }
   }
   Ket_occ_ph = arma::vec(occvec);
   Ket_unocc_ph = arma::vec(unoccvec);
}


//int TwoBodyChannel::GetLocalIndex(int p, int q) const { return KetMap[modelspace->GetKetIndex(p,q)];}; 
int TwoBodyChannel::GetLocalIndex(int p, int q) const
{
 if (p<=q)
   return KetMap[modelspace->GetKetIndex(p,q)];
 else
   return KetMap[modelspace->GetKetIndex(q,p)] + NumberKets;
} 

// get pointer to ket using local index
Ket & TwoBodyChannel::GetKet(int i) { return modelspace->GetKet(KetList[i]);}; 


//bool TwoBodyChannel::CheckChannel_ket(int p, int q) const
bool TwoBodyChannel::CheckChannel_ket(Orbit* op, Orbit* oq) const
{
   if ((op->index==oq->index) and (J%2 != 0)) return false; // Pauli principle
   if ((op->l + oq->l)%2 != parity) return false;
   if ((op->tz2 + oq->tz2) != 2*Tz) return false;
   if (op->j2 + oq->j2 < 2*J)       return false;
   if (abs(op->j2 - oq->j2) > 2*J)  return false;

   return true;
}

arma::uvec& TwoBodyChannel::GetKetIndex_pp() { return KetIndex_pp;};
arma::uvec& TwoBodyChannel::GetKetIndex_hh() { return KetIndex_hh;};
arma::uvec& TwoBodyChannel::GetKetIndex_ph() { return KetIndex_ph;};
arma::uvec& TwoBodyChannel::GetKetIndex_cc() { return KetIndex_cc;};
arma::uvec& TwoBodyChannel::GetKetIndex_vc() { return KetIndex_vc;};
arma::uvec& TwoBodyChannel::GetKetIndex_qc() { return KetIndex_qc;};
arma::uvec& TwoBodyChannel::GetKetIndex_vv() { return KetIndex_vv;};
arma::uvec& TwoBodyChannel::GetKetIndex_qv() { return KetIndex_qv;};
arma::uvec& TwoBodyChannel::GetKetIndex_qq() { return KetIndex_qq;};



arma::uvec TwoBodyChannel::GetKetIndexFromList(vector<index_t>& vec_in)
{
   vector<index_t> index_list (min(vec_in.size(),KetList.size()));
   auto it = set_intersection(KetList.begin(),KetList.end(),vec_in.begin(),vec_in.end(),index_list.begin());
   index_list.resize(it-index_list.begin());
   for (auto& x : index_list)
   {
     x = KetMap[x];
   }
   return arma::uvec(index_list);
}

//************************************************************************
//************************************************************************
//************************************************************************

TwoBodyChannel_CC::~TwoBodyChannel_CC()
{
//   cout << "In TwoBodyChannel_CC destructor" << endl;
}

TwoBodyChannel_CC::TwoBodyChannel_CC()
{}

TwoBodyChannel_CC::TwoBodyChannel_CC(int j, int p, int t, ModelSpace *ms)
{
  Initialize(ms->GetTwoBodyChannelIndex(j,p,t), ms);
}

TwoBodyChannel_CC::TwoBodyChannel_CC(int N, ModelSpace *ms)
{
   Initialize(N,ms);
}

// Check if orbits pq participate in this cross-coupled two-body channel
// Difference from regular channels:
// no Pauli rule, <pp||nn> is allowed. But |Tz| is still conserved,
// i.e. <pp||pn> is not allowed. So we use |Tz| rather than Tz,
// and don't use Tz=-1.
bool TwoBodyChannel_CC::CheckChannel_ket(Orbit* op, Orbit* oq) const
{
   if ((op->l + oq->l)%2 != parity)    return false;
   if (abs(op->tz2 + oq->tz2) != 2*Tz) return false;
   if (op->j2 + oq->j2 < 2*J)          return false;
   if (abs(op->j2 - oq->j2) > 2*J)     return false;

   return true;
}


//************************************************************************
//************************************************************************

// Static members
unordered_map<long long unsigned int, long double> ModelSpace::radList;
unordered_map<long long unsigned int, long double> ModelSpace::OsToHydroCoeffList;
unordered_map<unsigned long int,double> ModelSpace::SixJList;
unordered_map<unsigned long long int,double> ModelSpace::NineJList;
unordered_map<unsigned long long int,double> ModelSpace::MoshList;
map<string,vector<string>> ModelSpace::ValenceSpaces  {
{ "s-shell"  ,         {"vacuum", "p0s1","n0s1"}},
{ "p-shell"  ,         {"He4", "p0p3","n0p3","p0p1","n0p1"}},
{ "sp-shell"  ,        {"vacuum", "p0s1","n0s1","p0p3","n0p3","p0p1","n0p1"}},
{ "sd-shell"  ,        {"O16", "p0d5","n0d5","p0d3","n0d3","p1s1","n1s1"}},
{ "psd-shell"  ,       {"He4", "p0p3","n0p3","p0p1","n0p1","p0d5","n0d5","p0d3","n0d3","p1s1","n1s1"}},
{ "psdNR-shell"  ,     {"He10","p0p3","p0p1","n0d5","n0d3","n1s1"}}, // protons in p shell, neutrons in sd shell (NR is for neutron-rich)
{ "psdPR-shell"  ,     {"O10","n0p3","n0p1","p0d5","p0d3","p1s1"}}, // neutrons in p shell, protnons in sd shell (PR is for proton-rich)
{ "fp-shell"  ,        {"Ca40","p0f7","n0f7","p0f5","n0f5","p1p3","n1p3","p1p1","n1p1"}},
{ "sdfp-shell"  ,      {"O16", "p0d5","n0d5","p0d3","n0d3","p1s1","n1s1","p0f7","n0f7","p0f5","n0f5","p1p3","n1p3","p1p1","n1p1"}},
{ "sdfpNR-shell"  ,    {"O28", "p0d5","p0d3","p1s1","n0f7","n0f5","n1p3","n1p1"}}, // protons in sd shell, neutrons in fp shell, a la SDPFU from Nowacki/Poves (NR is for neutron-rich)
{ "sdfpPR-shell"  ,    {"Ca28", "n0d5","n0d3","n1s1","p0f7","p0f5","p1p3","p1p1"}}, // neutrons in sd shell, protons in fp shell, (PR is for proton-rich)
{ "fpg9-shell"  ,      {"Ca40","p0f7","n0f7","p0f5","n0f5","p1p3","n1p3","p1p1","n1p1","p0g9","n0g9"}},
{ "fpg9NR-shell"  ,    {"Ca40","p0f7","n0f7","p0f5","n0f5","p1p3","n1p3","p1p1","n1p1","n0g9"}}, // just add g9/2 for neutrons
{ "fpgdsNR-shell"  ,   {"Ca60","p0f7","p0f5","p1p3","p1p1","n0g9","n0g7","n1d5","n1d3","n2s1"}}, // protons in the fp shell, neutrons in the gds shell
//{ "fpg9NR-shell"  ,    {"Ca52","p0f7","p0f5","p1p3","p1p1","n0f5","n1p1","n0g9"}}, // protons in the fp shell, neutrons in upper fp + g9/2
{ "sd3f7p3-shell"  ,   {"Si28","p0d3","n0d3","p1s1","n1s1","p0f7","n0f7","p1p3","n1p3"}},
{ "gds-shell" ,        {"Zr80","p0g9","n0g9","p0g7","n0g7","p1d5","n1d5","p1d3","n1d3","p2s1","n2s1"}}, // This is a big valence space, more than a few particles will be a serious shell model diagonalization
};



ModelSpace::~ModelSpace()
{
//  cout << "In ModelSpace destructor. emax = " << Emax << endl;
}

ModelSpace::ModelSpace()
:  Emax(0), E2max(0), E3max(0), Lmax(0), Lmax2(0), Lmax3(0), OneBodyJmax(0), TwoBodyJmax(0), ThreeBodyJmax(0), norbits(0),
  hbar_omega(20), target_mass(16), SystemType("nuclear")
{
  cout << "In default constructor" << endl;
}


ModelSpace::ModelSpace(const ModelSpace& ms)
 :
   holes( ms.holes), particles( ms.particles),
   core(ms.core), valence(ms.valence), qspace( ms.qspace), 
   proton_orbits( ms.proton_orbits),neutron_orbits( ms.neutron_orbits),
   KetIndex_pp( ms.KetIndex_pp), KetIndex_ph( ms.KetIndex_ph), KetIndex_hh( ms.KetIndex_hh),
   KetIndex_cc( ms.KetIndex_cc),
   KetIndex_vc( ms.KetIndex_vc),
   KetIndex_qc( ms.KetIndex_qc),
   KetIndex_vv( ms.KetIndex_vv),
   KetIndex_qv( ms.KetIndex_qv),
   KetIndex_qq( ms.KetIndex_qq),
   Ket_occ_hh( ms.Ket_occ_hh),
   Ket_unocc_hh( ms.Ket_unocc_hh),
   Emax(ms.Emax), E2max(ms.E2max), E3max(ms.E3max), Lmax(ms.Lmax), Lmax2(ms.Lmax2), Lmax3(ms.Lmax3),
   OneBodyJmax(ms.OneBodyJmax), TwoBodyJmax(ms.TwoBodyJmax), ThreeBodyJmax(ms.ThreeBodyJmax),
   OneBodyChannels(ms.OneBodyChannels),
   SortedTwoBodyChannels(ms.SortedTwoBodyChannels),
   SortedTwoBodyChannels_CC(ms.SortedTwoBodyChannels_CC),
   norbits(ms.norbits), hbar_omega(ms.hbar_omega),
   target_mass(ms.target_mass), target_Z(ms.target_Z), Aref(ms.Aref), Zref(ms.Zref),
   nTwoBodyChannels(ms.nTwoBodyChannels),
   Orbits(ms.Orbits), Kets(ms.Kets),
   TwoBodyChannels(ms.TwoBodyChannels), TwoBodyChannels_CC(ms.TwoBodyChannels_CC),
   SystemType(ms.SystemType),
   systemBasis(ms.systemBasis)
{
   for (TwoBodyChannel& tbc : TwoBodyChannels)   tbc.modelspace = this;
   for (TwoBodyChannel_CC& tbc_cc : TwoBodyChannels_CC)   tbc_cc.modelspace = this;
}

ModelSpace::ModelSpace(ModelSpace&& ms)
 :
   holes( move(ms.holes)), particles( move(ms.particles)),
   core(move(ms.core)), valence(move(ms.valence)),  qspace( move(ms.qspace)),  
   proton_orbits( move(ms.proton_orbits)),
   neutron_orbits( move(ms.neutron_orbits)),
   KetIndex_pp( move(ms.KetIndex_pp)), KetIndex_ph( move(ms.KetIndex_ph)), KetIndex_hh( move(ms.KetIndex_hh)),
   KetIndex_cc( ms.KetIndex_cc),
   KetIndex_vc( ms.KetIndex_vc),
   KetIndex_qc( ms.KetIndex_qc),
   KetIndex_vv( ms.KetIndex_vv),
   KetIndex_qv( ms.KetIndex_qv),
   KetIndex_qq( ms.KetIndex_qq),
   Ket_occ_hh( ms.Ket_occ_hh),
   Ket_unocc_hh( ms.Ket_unocc_hh),
   Emax(ms.Emax), E2max(ms.E2max), E3max(ms.E3max), Lmax(ms.Lmax), Lmax2(ms.Lmax2), Lmax3(ms.Lmax3),
   OneBodyJmax(ms.OneBodyJmax), TwoBodyJmax(ms.TwoBodyJmax), ThreeBodyJmax(ms.ThreeBodyJmax),
   OneBodyChannels(move(ms.OneBodyChannels)),
   SortedTwoBodyChannels(move(ms.SortedTwoBodyChannels)),
   SortedTwoBodyChannels_CC(move(ms.SortedTwoBodyChannels_CC)),
   norbits(ms.norbits), hbar_omega(ms.hbar_omega),
   target_mass(ms.target_mass), target_Z(ms.target_Z), Aref(ms.Aref), Zref(ms.Zref),
   nTwoBodyChannels(ms.nTwoBodyChannels),
   Orbits(move(ms.Orbits)), Kets(move(ms.Kets)),
   TwoBodyChannels(move(ms.TwoBodyChannels)), TwoBodyChannels_CC(move(ms.TwoBodyChannels_CC)),
   SystemType(move(ms.SystemType)),
   systemBasis(move(ms.systemBasis))
{
   for (TwoBodyChannel& tbc : TwoBodyChannels)   tbc.modelspace = this;
   for (TwoBodyChannel_CC& tbc_cc : TwoBodyChannels_CC)   tbc_cc.modelspace = this;
   for (TwoBodyChannel& tbc : ms.TwoBodyChannels)   tbc.modelspace = NULL;
   for (TwoBodyChannel_CC& tbc_cc : ms.TwoBodyChannels_CC)   tbc_cc.modelspace = NULL;
}


// orbit string representation is e.g. p0f7
// Assumes that the core is hole states that aren't in the valence space.
ModelSpace::ModelSpace(int emax, vector<string> hole_list, vector<string> valence_list, int Lmax, string SystemType, string systemBasis)
:  Emax(emax), E2max(2*emax), E3max(3*emax), Lmax(Lmax), Lmax2(emax), Lmax3(emax), OneBodyJmax(0), TwoBodyJmax(0), ThreeBodyJmax(0), norbits(0), hbar_omega(20), target_mass(16)
{
   Init(emax, hole_list, hole_list, valence_list,Lmax, SystemType, systemBasis); 
}

// If we don't want the reference to be the core
ModelSpace::ModelSpace(int emax, vector<string> hole_list, vector<string> core_list, vector<string> valence_list, int Lmax, string SystemType, string systemBasis)
: Emax(emax), E2max(2*emax), E3max(3*emax), Lmax(Lmax), Lmax2(emax), Lmax3(emax), OneBodyJmax(0), TwoBodyJmax(0), ThreeBodyJmax(0), norbits(0), hbar_omega(20), target_mass(16)
{
   Init(emax, hole_list, core_list, valence_list,Lmax, SystemType, systemBasis); 
}

// Most conventient interface
ModelSpace::ModelSpace(int emax, string reference, string valence, int Lmax, string SystemType, string systemBasis)
: Emax(emax), E2max(2*emax), E3max(3*emax), Lmax(Lmax), Lmax2(emax), Lmax3(emax), OneBodyJmax(0), TwoBodyJmax(0), ThreeBodyJmax(0),hbar_omega(20)
{
  Init(emax,reference,valence,Lmax, SystemType, systemBasis); // Need to ensure all of these are init'd before passed.
}

ModelSpace::ModelSpace(int emax, string valence, int Lmax, string SystemType, string systemBasis)
: Emax(emax), E2max(2*emax), E3max(3*emax), Lmax(Lmax), Lmax2(emax), Lmax3(emax), OneBodyJmax(0), TwoBodyJmax(0), ThreeBodyJmax(0),hbar_omega(20),SystemType("nuclear"), systemBasis("harmonic")
{
  auto itval = ValenceSpaces.find(valence);
  if ( itval != ValenceSpaces.end() ) // we've got a valence space
     Init(emax,itval->second[0],valence,Lmax, systemBasis);
  else  // no valence space. we've got a single-reference.
     Init(emax,valence,valence,Lmax, SystemType, systemBasis);
}



// Specify the reference and either the core or valence
// This is the most convenient interface
void ModelSpace::Init(int emax, string reference, string valence, int Lmax, string SystemType, string systemBasis)
{
//  int Aref,Zref;
  cout << "Building ModelSpace..." << endl;
  GetAZfromString(reference,Aref,Zref);
  if (SystemType == "nuclear"){
	map<index_t,double> hole_list = GetOrbitsAZ(Aref,Zref);
	Init(emax,hole_list,valence,Lmax, SystemType, systemBasis);}
  else if (SystemType == "atomic"){
	cout << "Atomic basis selected, getting hole_list..." << endl;
	map<index_t,double> hole_list = GetOrbitsE(Zref); // This was an attempt to redo the way modelspaces are built
	//map<index_t,double> hole_list = GetOrbitsAZ(Aref,Zref);
	Init(emax,hole_list,valence,Lmax, SystemType, systemBasis);} // Need to ensure all of these are init'd before passed.
}

void ModelSpace::Init(int emax, map<index_t,double> hole_list, string valence, int Lmax, string SystemType, string systemBasis)
{
  int Ac,Zc;
  vector<index_t> valence_list, core_list;
//  map<index_t,double> core_map;

  if (valence == "0hw-shell")
  {
    Get0hwSpace(Aref,Zref,core_list,valence_list);
  }
  else
  {
    auto itval = ValenceSpaces.find(valence);
  
    if ( itval != ValenceSpaces.end() ) // we've got a valence space
    {
       string core_string = itval->second[0];
       GetAZfromString(core_string,Ac,Zc);
       valence_list = String2Index(vector<string>(itval->second.begin()+1,itval->second.end()));
    }
    else  // no valence space. we've got a single-reference.
    {
       GetAZfromString(valence,Ac,Zc);
    }
  
    //core_map = GetOrbitsAZ(Ac,Zc);
    //for (auto& it_core : core_map) core_list.push_back(it_core.first);
    cout << "Getting core_list..." << endl;
    if (SystemType == "atomic") for (auto& it_core : GetOrbitsE(Zc) ) core_list.push_back(it_core.first);
    if (SystemType == "nuclear") for (auto& it_core : GetOrbitsAZ(Ac,Zc) ) core_list.push_back(it_core.first);
  }

  target_mass = Aref;
  target_Z = Zref;
  Init(emax, hole_list,core_list, valence_list, Lmax, SystemType, systemBasis); // Need to ensure all of these are init'd before passed.
  
}


// Specify the model space with strings of orbit lists.
// Less convenient, but more flexible
void ModelSpace::Init(int emax, vector<string> hole_list, vector<string> core_list, vector<string> valence_list, int Lmax, string SystemType, string systemBasis)
{
   cout << "Creating a model space with Emax = " << Emax << "  and hole orbits [";
   for (auto& h : hole_list)  cout << h << " ";
   cout << "]   and core orbits [";
   for (auto& c : core_list)    cout << c << " ";
   cout << "]   and valence orbits [";
   for (auto& v : valence_list)   cout << v << " ";
   cout << "]" << endl;
   map<index_t,double> hole_map;
   for (auto& h : String2Index(hole_list)) hole_map[h] = 1.0;
  Init(emax, hole_map, String2Index(core_list), String2Index(valence_list), Lmax, SystemType, systemBasis);
}


void ModelSpace::Init_occ_from_file(int emax, string valence, string occ_file, int Lmax, string SystemType, string systemBasis)
{
  index_t orb;
  double occ;
  map<index_t,double> hole_list;

  ifstream infile(occ_file);
  if (!infile.good())
  {
    cout << endl << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
    cout << "Trouble reading file: " << occ_file << endl;
    cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl << endl;
  }

  while( infile >> orb >> occ )
  {
    if ( hole_list.find(orb) != hole_list.end() and  abs( hole_list[orb] -occ) > -1e-6) // the minus sign is for a test. Change it back.
    {
        cout << "Warning: in file " << occ_file << ", redefinition of occupation of orbit "
             << orb << "  " << hole_list[orb] << " => " << occ << endl;
    }
    cout << "from occ file: " << endl;
    hole_list[orb] = occ;
    cout << orb << " " << occ << endl;
  }

  Init(emax, hole_list, valence, Lmax, SystemType, systemBasis);
}


// This is the Init which should inevitably be called
void ModelSpace::Init(int emax, map<index_t,double> hole_list, vector<index_t> core_list, vector<index_t> valence_list, int Lmax, string SystemType, string systemBasis)
{
   cout << "SystemBasis=" << systemBasis << endl;
   ClearVectors();
   if (Lmax < 0) Lmax = emax;
   emax = Emax;
   cout << "core list: ";
   for (auto& c : core_list) cout << c << " ";
   cout << endl;
   cout << "valence list: ";
   for (auto& v : valence_list) cout << v << " ";
   cout << endl;
   cout << "hole list: ";
   for (auto& h : hole_list) cout << h.first << " ( " << h.second << " ) ";
   cout << endl;

   // Make sure no orbits are both core and valence
   for (auto& c : core_list)
   {
     if ( find(valence_list.begin(), valence_list.end(), c) != valence_list.end() )
       cout << "!!!!!!!!!!!!! ModelSpace::Init : Conflicting definition. Orbit " << c << " is in core and valence spaces." << endl;
   }

   norbits = (Emax+1)*(Emax+2); // Need to take into account effect of Lmax
   Orbits.resize(norbits);
   int real_norbits = 0;
   int count = 0;
   if (systemBasis == "harmonic"){
       for (int N=0; N<=Emax; ++N)
       {
         //min(N,Lmax)
         //cout << "Lmax=" << Lmax << " N=" << N << endl;
	 int lmax = 0;
	 int lmin = 0;
	 int ldiff = 0;
	 if (SystemType == "atomic")
	   {
		lmax = N;
		if (N%2==0) lmin=0;
		else lmin=1;
		ldiff = 2;
	   } else {
		lmax = 0;
		lmin = N;
		ldiff = -2;
	   }
         for (int l=lmin; l<=lmax; l+=ldiff)
         {
           if (l>Lmax) continue;
           int n = (N-l)/2;
	   int jmax = 0;
	   int jmin = 0;
	   int jdiff = 0;
	   if (SystemType == "atomic")
	   {
		jmax = 2*l+1;
		jmin = abs(2*l-1);
		jdiff = 2;
	   } else {
		jmax = abs(2*l-1);
		jmin = 2*l+1;
		jdiff = -2;
	   }
           for (int j2=jmin; j2<=jmax and j2>0; j2+=jdiff)
           {
             for (int tz : {-1, 1} )
             {
                double occ = 0;
                int cvq = 2;
                int indx = Index1(n,l,j2,tz);
    	    	if (SystemType == "atomic" and tz < 0){ // Checks twice, this looks like garbage
    		    indexMap[indx] = count; // Map atomic orbits as well as nuclear
    		    indx = indexMap[indx]; 
    		    //indexMap[count] = count;
    		    //indx = count;
		    count++;
		    //indx /= 2;
	        }
                if (hole_list.find(indx) != hole_list.end()) occ = hole_list[indx];
                if ( find(core_list.begin(), core_list.end(), indx) != core_list.end() ) cvq=0; // core orbit
                if ( find(valence_list.begin(), valence_list.end(), indx) != valence_list.end() ) cvq=1; // valence orbit
                cout << "Orbit with n=" << n << " l=" << l << " j2=" << j2 << " tz=" << tz << " occ=" << occ << " cvq=" << cvq << " at indx=" << indx << endl;
	        if (SystemType == "atomic" and tz < 0) { // Only allow isospin -1/2.  this simulates only protons being added
		    AddOrbit(n,l,j2,tz,occ,cvq,indx);
		    cout << "Added orbit." << endl;
		    real_norbits++;
	        } else if (SystemType == "nuclear" ) {
		    AddOrbit(n,l,j2,tz,occ,cvq);
		    real_norbits++;
	        } else {
		    //offset++; // In case we don't add an orbit, iterate a bit further to try to get all the orbits in.
	        }
             }
           }
         }
       }
   } else if (systemBasis == "hydrogen") {
	for (int n=1; n<=Emax; n++) 
	{
	    for (int l=0; l<=Lmax and l<n; l++)
	    {
		for (int j2=abs(2*l-1); j2 <= 2*l+1; j2+=2)
		{
		    int tz = -1;
		    double occ = 0;
		    int cvq = 2;
		    int indx = Index_atomic(n,l,j2); //Index1(n,l,j2,tz);
		    indexMap[indx] = count;
		    indx = indexMap[indx];
		    count++;
		    if (hole_list.find(indx) != hole_list.end()) occ = hole_list[indx];
		    if ( find(core_list.begin(), core_list.end(), indx) != core_list.end() ) cvq=0; // core orbit
		    if ( find(valence_list.begin(), valence_list.end(), indx) != valence_list.end() ) cvq=1; // valence orbit
		    cout << "Hydrogen Orbit with n=" << n << " l=" << l << " j2=" << j2 << " tz=" << tz << " occ=" << occ << " cvq=" << cvq << " at indx=" << indx << endl;
		    AddOrbit(n,l,j2,tz,occ,cvq,indx);
		    real_norbits++;
		}
	    }
	}
   }
	
   norbits = real_norbits;
   cout << "Reset norbits=" << norbits << endl;
   Orbits.resize(norbits);
   cout << "Resized Orbits; Orbits.size()=" << Orbits.size() << endl;
   Aref = 0;
   Zref = 0;
   for (auto& h : holes)
   {
     Orbit& oh = GetOrbit(h);
     Aref += (oh.j2+1)*oh.occ;
     if (oh.tz2 < 0) Zref += (oh.j2+1)*oh.occ;
   }
   cout << "Setting up kets." << endl;
   if (SystemType != ""){
      SetupKets(SystemType);
   } else {
      SetupKets("nuclear");
   }
   cout << "Set up kets." << endl;
}



// Get vector of orbit indices from vector of strings
// e.g. "p0f7" gives the index of the proton 0f7/2 orbit.
vector<index_t> ModelSpace::String2Index( vector<string> vs )
{
  vector<index_t> vi;
  vector<char> l_list = {'s','p','d','f','g','h','i','j','k','l','m','n','o'};

  for ( auto& s : vs )
  {
    int n,l,j2,tz2;
    tz2 = s[0]=='p' ? -1 : 1;
    istringstream( s.substr(1,2) ) >> n;
    l = find(l_list.begin(),l_list.end(), s[2]) - l_list.begin();
    istringstream( s.substr(3,s.size()) ) >> j2;
    vi.push_back( Index1(n,l,j2,tz2) );
  }
  return vi;
}


string ModelSpace::Index2String( index_t ind)
{
  vector<char> l_list = {'s','p','d','f','g','h','i','j','k','l','m','n','o'};
  Orbit& oi = GetOrbit(ind);
  char c[10];
  char pn = oi.tz2 < 0 ? 'p' : 'n';
  char lstr = l_list[oi.l];
  sprintf(c, "%c%d%c%d",pn,oi.n,lstr,oi.j2);
  return string(c) ;
}


void ModelSpace::GetAZfromString(string str,int& A, int& Z) // TODO: accept different formats, e.g. 22Na vs Na22
{
  vector<string> periodic_table = {"n","H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar",
                        "K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr",
                        "Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe",
                        "Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf",
                        "Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb"};
  if (str == "vacuum") str="n0";
  int i=0;
  while (! isdigit(str[i])) i++;
  string elem = str.substr(0,i);
  stringstream( str.substr(i,str.size()-i)) >> A;
  auto it_elem = find(periodic_table.begin(),periodic_table.end(),elem);
  if (it_elem != periodic_table.end())
  {
    Z = it_elem - periodic_table.begin();
  }
  else
  {
    Z =-1;
   cout << "ModelSpace::GetAZfromString :  Trouble parsing " << str << endl;
  }
}

// Returns the mapped orbit back to the function in the event that the system is atomic.

Orbit& ModelSpace::GetOrbit(int i){
   //if(SystemType == "atomic"){
   //   return (Orbit&) Orbits[indexMap[i]];
   //} else {
      return (Orbit&) Orbits[i];
   //}
}

// Fill A orbits with Z protons and A-Z neutrons
// assuming a standard shell-model level ordering
map<index_t,double> ModelSpace::GetOrbitsAZ(int A, int Z)
{
  int zz = 0;
  int nn = 0; // unfortunate there are so many n's here...
  map<index_t,double> holesAZ;
  for (int N=0; N<=Emax; ++N)
  {
    for (int g=2*N+1;g>=-2*N;g-=4)
    {
      int j2 = abs(g);
      int l = g<0 ? (j2+1)/2 : (j2-1)/2;
      int n = (N-l)/2;

      if (zz < Z)
      {
        int dz = min(Z-zz,j2+1);
	int indx = Index1(n,l,j2,-1);
	//if (SystemType == "atomic"){ 
	//  indx = indexMap[indx];}
        holesAZ[indx] = dz/(j2+1.0);
        zz += dz;
      }
      if (nn < A-Z)
      {
        int dn = min(A-Z-nn,j2+1);
        holesAZ[Index1(n,l,j2,1)] = dn/(j2+1.0);
        nn += dn;
      }
      if (zz==Z and nn==A-Z)
      {
         return holesAZ; // We're all done here.
      }
    }
  }
  cout << "Trouble! Model space not big enough to fill A=" << A << " Z="<< Z << "  emax = " << Emax << endl;
  return holesAZ;

}

map<index_t,double> ModelSpace::GetOrbitsE(int Z)
{
    int z=0;
    map<index_t,double> holesE;
    cout << "In GetOrbitsE for Z=" << Z << endl;
   for (int N=1; N<=Emax; ++N)
   {
	//cout << "Gettin' ready!" << endl;
	for (int l=0; l <= N and l <= Lmax; l++)
	{
	    for (int j2=abs(2*l-1); j2 <= 2*l +1; j2+=2)
	    {
		int n = N; //(N-l)/2;
		int indx = Index_atomic(n,l,j2);
		cout << "N=" << N << " j2=" << j2 << " l=" << l << " n=" << n << " Index1=" << indx << endl;
		if (z < Z)
		{
		    int dz = min(Z-z, j2+1);
		    cout << "dz=" << dz << " z=" << z << " dz/(j2+1.0)=" << dz/(j2+1.0) << endl;
		    //indexMap[indx] = indexMap.size()-1;
		    //holesE[indexMap[Index1(n,l,j2,-1)]] = dz/(j2+1.0);
		    holesE[indx] = dz/(j2+1.0); // indx/2 ?
		    z += dz;
		}
	    }
	}
	if (z == Z) return holesE; // We're all done here.
    }
    cout << "Didn't set ModelSpace big enough to fill Z=" << Z << " with emax = " << Emax << endl;
    return holesE;
}


/// Find the valence space of one single major oscillator shell each for protons and neutrons (not necessarily
/// the same shell for both) which contains the naive shell-model ground state of the reference.
/// For example, if we want to treat C20, with 6 protons and 14 neutrons, we take the 0p shell for protons
/// and 1s0d shell for neutrons.
void ModelSpace::Get0hwSpace(int Aref, int Zref, vector<index_t>& core_list, vector<index_t>& valence_list)
{
  int Nref = Aref-Zref;
  int OSC_protons=0,OSC_neutrons=0;
  while ( (OSC_protons +1)*(OSC_protons +2)*(OSC_protons +3)/3 <= Zref ) OSC_protons++;
  while ( (OSC_neutrons+1)*(OSC_neutrons+2)*(OSC_neutrons+3)/3 <= Nref ) OSC_neutrons++;

  int Zcore = (OSC_protons )*(OSC_protons +1)*(OSC_protons +2)/3;
  int Ncore = (OSC_neutrons)*(OSC_neutrons+1)*(OSC_neutrons+2)/3;

  for (auto& it_core : GetOrbitsAZ(Zcore+Ncore,Zcore)) core_list.push_back(it_core.first);

  for (int L=OSC_protons; L>=0; L-=2)
  {
    for (int j2=2*L+1;j2>max(2*L-2,0);j2-=2)
    {
      valence_list.push_back( GetOrbitIndex( (OSC_protons-L)/2, L, j2, -1) );
    }
  }
  for (int L=OSC_neutrons; L>=0; L-=2)
  {
    for (int j2=2*L+1;j2>max(2*L-2,0);j2-=2)
    {
      valence_list.push_back( GetOrbitIndex( (OSC_neutrons-L)/2, L, j2, 1) );
    }
  }

}



void ModelSpace::SetReference(vector<index_t> new_reference)
{
  vector<index_t> c = core;
  vector<index_t> v = valence;
  map<index_t,double> h;
  for (auto r : new_reference) h[r] = 1.0;
  ClearVectors();
  Init(Emax, h,c,v);
}

void ModelSpace::SetReference(map<index_t,double> new_reference)
{
  vector<index_t> c = core;
  vector<index_t> v = valence;
  ClearVectors();
  Init(Emax, new_reference,c,v);
}

void ModelSpace::SetReference(string new_reference)
{
  vector<index_t> c = core;
  vector<index_t> v = valence;
  ClearVectors();
  GetAZfromString(new_reference,Aref,Zref);
  map<index_t,double> h = GetOrbitsAZ(Aref,Zref);
  Init(Emax, h,c,v);
}

ModelSpace ModelSpace::operator=(const ModelSpace& ms)
{
   holes =  ms.holes;
   particles =  ms.particles;
   valence = ms.valence;
   qspace =  ms.qspace;
   core = ms.core;
   proton_orbits =  ms.proton_orbits;
   neutron_orbits =  ms.neutron_orbits;
   KetIndex_pp =  ms.KetIndex_pp;
   KetIndex_ph =  ms.KetIndex_ph;
   KetIndex_hh =  ms.KetIndex_hh;
   KetIndex_cc =  ms.KetIndex_cc;
   KetIndex_vc =  ms.KetIndex_vc;
   KetIndex_qc =  ms.KetIndex_qc;
   KetIndex_vv =  ms.KetIndex_vv;
   KetIndex_qv =  ms.KetIndex_qv;
   KetIndex_qq =  ms.KetIndex_qq;
   Ket_occ_hh  =  ms.Ket_occ_hh;
   Ket_unocc_hh  =  ms.Ket_unocc_hh;
   Emax = ms.Emax;
   E2max = ms.E2max;
   E3max = ms.E3max;
   Lmax2 = ms.Lmax2;
   Lmax3 = ms.Lmax3;
   OneBodyJmax = ms.OneBodyJmax;
   TwoBodyJmax = ms.TwoBodyJmax;
   ThreeBodyJmax = ms.ThreeBodyJmax;
   OneBodyChannels = ms.OneBodyChannels;
   SortedTwoBodyChannels = ms.SortedTwoBodyChannels;
   SortedTwoBodyChannels_CC = ms.SortedTwoBodyChannels_CC;
   norbits = ms.norbits;
   hbar_omega = ms.hbar_omega;
   target_mass = ms.target_mass;
   target_mass = ms.target_Z;
   Aref = ms.Aref;
   Zref = ms.Zref;
   Orbits = ms.Orbits;
   Kets = ms.Kets;
   TwoBodyChannels = ms.TwoBodyChannels;
   TwoBodyChannels_CC = ms.TwoBodyChannels_CC;
   for (TwoBodyChannel& tbc : TwoBodyChannels)   tbc.modelspace = this;
   for (TwoBodyChannel_CC& tbc_cc : TwoBodyChannels_CC)   tbc_cc.modelspace = this;

//   cout << "In copy assignment for ModelSpace" << endl;
   return ModelSpace(*this);
}



ModelSpace ModelSpace::operator=(ModelSpace&& ms)
{
   holes =  move(ms.holes);
   particles =  move(ms.particles);
   valence = move(ms.valence);
   qspace =  move(ms.qspace);
   core = move(ms.core);
   proton_orbits =  move(ms.proton_orbits);
   neutron_orbits =  move(ms.neutron_orbits);
   KetIndex_pp =  move(ms.KetIndex_pp);
   KetIndex_ph =  move(ms.KetIndex_ph);
   KetIndex_hh =  move(ms.KetIndex_hh);
   KetIndex_cc =  move(ms.KetIndex_cc);
   KetIndex_vc =  move(ms.KetIndex_vc);
   KetIndex_qc =  move(ms.KetIndex_qc);
   KetIndex_vv =  move(ms.KetIndex_vv);
   KetIndex_qv =  move(ms.KetIndex_qv);
   KetIndex_qq =  move(ms.KetIndex_qq);
   Ket_unocc_hh =  move(ms.Ket_unocc_hh);
   Ket_occ_hh =  move(ms.Ket_occ_hh);
   Emax = move(ms.Emax);
   E2max = move(ms.E2max);
   E3max = move(ms.E3max);
   Lmax2 = move(ms.Lmax2);
   Lmax3 = move(ms.Lmax3);
   OneBodyJmax = move(ms.OneBodyJmax);
   TwoBodyJmax = move(ms.TwoBodyJmax);
   ThreeBodyJmax = move(ms.ThreeBodyJmax);
   OneBodyChannels = move(ms.OneBodyChannels);
   SortedTwoBodyChannels = move(ms.SortedTwoBodyChannels);
   SortedTwoBodyChannels_CC = move(ms.SortedTwoBodyChannels_CC);
   norbits = move(ms.norbits);
   hbar_omega = move(ms.hbar_omega);
   target_mass = move(ms.target_mass);
   target_Z = move(ms.target_Z);
   Aref = move(ms.Aref);
   Zref = move(ms.Zref);
   Orbits = move(ms.Orbits);
   Kets = move(ms.Kets);
   TwoBodyChannels = move(ms.TwoBodyChannels);
   TwoBodyChannels_CC = move(ms.TwoBodyChannels_CC);
   for (TwoBodyChannel& tbc : TwoBodyChannels)   tbc.modelspace = this;
   for (TwoBodyChannel_CC& tbc_cc : TwoBodyChannels_CC)   tbc_cc.modelspace = this;
   for (TwoBodyChannel& tbc : ms.TwoBodyChannels)   tbc.modelspace = NULL;
   for (TwoBodyChannel_CC& tbc_cc : ms.TwoBodyChannels_CC)   tbc_cc.modelspace = NULL;
   return ModelSpace(*this);
}



void ModelSpace::AddOrbit(Orbit orb)
{
  AddOrbit(orb.n, orb.l, orb.j2, orb.tz2, orb.occ, orb.cvq);
}

void ModelSpace::AddOrbit(int n, int l, int j2, int tz2, double occ, int cvq, int index)
{
   int ind = index;
   if (index == -1)
	index = Index1(n, l, j2, tz2);
   else
	ind = index;
   //index_t ind = Index1(n, l, j2, tz2)
   //ind = indexMap[ind];
   Orbits[index] = Orbit(n,l,j2,tz2,occ,cvq,index);
   if (j2 > OneBodyJmax)
   {
      OneBodyJmax = j2;
      TwoBodyJmax = OneBodyJmax;
      ThreeBodyJmax = OneBodyJmax*3-1;
      nTwoBodyChannels = 2*3*(TwoBodyJmax+1);
   }


   if ( occ < OCC_CUT) particles.push_back(ind);
   else holes.push_back(ind);
   if (cvq == 0) core.push_back(ind);
   if (cvq == 1) valence.push_back(ind);
   if (cvq == 2) qspace.push_back(ind);
   if (tz2 < 0 ) proton_orbits.push_back(ind);
   if (tz2 > 0 ) neutron_orbits.push_back(ind);

   OneBodyChannels[{l, j2, tz2}].push_back(ind);
}



int ModelSpace::GetOrbitIndex(string orb)
{
  vector<char> l_list = {'s','p','d','f','g','h','i','j','k','l','m','n','o'};
  int n=-1,l=-1,j2=-1;
  int tz2 = orb[0]=='p' ? -1 : 1;
  stringstream(orb.substr(1,1)) >> n;
  auto it_l = find(l_list.begin(), l_list.end(), orb[2]);
  if ( it_l != l_list.end() )
    l = it_l - l_list.begin();
  else
    cout << "Bad orbit label " << orb << endl;
  stringstream(orb.substr(3)) >> j2;
  return Index1(n,l,j2,tz2);
}

int ModelSpace::GetTwoBodyChannelIndex(int j, int p, int t)
{
   return (t+1)*2*(TwoBodyJmax+1) + p*(TwoBodyJmax+1) + j;
}




void ModelSpace::SetupKets(string Sys)
{
   cout << "Entering SetupKets()" << endl;
   int index = 0;
   //if (SystemType == "nuclear")
   //{
      Kets.resize(Index2(norbits-1,norbits-1)+1);
   //} else {
   //   Kets.resize(0);
   //}
   //Kets.resize(0.5*norbits*norbits + 0.5*norbits);
   //int count = 0;
   for (int p=0;p<norbits;p++)
   {
     for (int q=p;q<norbits;q++)
     {
	if (Sys == "atomic")
	{
	   index = Index2(p,q);
	   //cout << "Grabbing ket with p=" << p << " q=" << q << " at indexMap[p]= " << indexMap[p] << " and indexMap[q]=" << indexMap[q] << "and setting to index=" << index << endl;
	   //index = Kets.size();
	   //Kets.emplace_back(Ket(GetOrbit(p),GetOrbit(q)));
	   Orbit& o1 = GetOrbit(p);
	   Orbit& o2 = GetOrbit(q);
	   //index = Index2(o1.index,o2.index);
	   Kets[index] = Ket(o1,o2);
	   //Kets[index] = Ket(GetOrbit(indexMap[p]),GetOrbit(indexMap[q]));
	   //count++;
	} else
	{
	   index = Index2(p,q);
	   Kets[index] = Ket(GetOrbit(p),GetOrbit(q));
	}
        //cout << "index=" << index << " p=" << p << " q=" << q << endl;
        //Orbit& orbp = GetOrbit(p);
	//cout << "orb(" << p << ") n=" << orbp.n << " l=" << orbp.l << " j2=" << orbp.j2 << " tz2=" << orbp.tz2 << endl;
	//Orbit& orbq = GetOrbit(q);
	//cout << "orb(" << q << ") n=" << orbq.n << " l=" << orbq.l << " j2=" << orbq.j2 << " tz2=" << orbq.tz2 << endl;
        
     }
   }
  
  //cout << "Set up Kets[], moving to Ket& ket; Kets[].size()=" << Kets.size() << endl;
  for (index_t index=0;index<Kets.size();++index)
  {
    //cout << "index=" << index << endl;
    Ket& ket = Kets[index];
    //cout << "Got the ket, checking parity and Tz." << endl;
    int Tz = (ket.op->tz2 + ket.oq->tz2)/2;
    int parity = (ket.op->l + ket.oq->l)%2;
    //cout << "ket.op->l=" << ket.op->l << " ket.oq->l=" << ket.oq->l << endl;
    //cout << "About to add MonopoleKet with Tz=" << Tz << " parity=" << parity << " at index=" << index << endl;
    MonopoleKets[Tz+1][parity][index] = MonopoleKets[Tz+1][parity].size()-1;
    //cout << "Added MonopoleKet." << endl;
    double occp = ket.op->occ;
    double occq = ket.oq->occ;
    int cvq_p = ket.op->cvq;
    int cvq_q = ket.oq->cvq;
    if (cvq_p+cvq_q==0)      KetIndex_cc.push_back(index); // 00
    if (cvq_p+cvq_q==1)      KetIndex_vc.push_back(index); // 01
    if (abs(cvq_p-cvq_q)==2) KetIndex_qc.push_back(index); // 02
    if (cvq_p*cvq_q==1)      KetIndex_vv.push_back(index); // 11
    if (cvq_p+cvq_q==3)      KetIndex_qv.push_back(index); // 12
    if (cvq_p+cvq_q==4)      KetIndex_qq.push_back(index); // 22
    if (occp<OCC_CUT and occq<OCC_CUT) KetIndex_pp.push_back(index);
//    if (occp>OCC_CUT or occq>OCC_CUT)
    if ( (occp>OCC_CUT) xor (occq>OCC_CUT) )
    {
       KetIndex_ph.push_back(index);
       Ket_occ_ph.push_back(occp*occq);
       Ket_unocc_ph.push_back((1-occp)*(1-occq));
    }
    if (occp>OCC_CUT and occq>OCC_CUT)
    {
       KetIndex_hh.push_back(index);
       Ket_occ_hh.push_back(occp*occq);
       Ket_unocc_hh.push_back((1-occp)*(1-occq));
    }
    //cout << "About to loop again." << endl;
   }
   //cout << "Got past Ket&; resizing TB." << endl;
   SortedTwoBodyChannels.resize(nTwoBodyChannels);
   SortedTwoBodyChannels_CC.resize(nTwoBodyChannels);
   //cout << "Resized TB; sorting TB., nTwoBodyChannels=" << nTwoBodyChannels << endl;
   for (int ch=0;ch<nTwoBodyChannels;++ch)
   {
      TwoBodyChannels.push_back(move(TwoBodyChannel(ch,this)));
      TwoBodyChannels_CC.push_back(move(TwoBodyChannel_CC(ch,this)));
      SortedTwoBodyChannels[ch] = ch;
      SortedTwoBodyChannels_CC[ch] = ch;
      //cout << "ch=" << ch << " nkets=" << TwoBodyChannels[ch].GetNumberKets() << endl;
      //cout << "CC_ch=" << ch << " CC_nkets=" << TwoBodyChannels_CC[ch].GetNumberKets() << endl;
   }
   //cout << "Initialized dem channels." << endl;
   // Sort the two body channels in descending order of matrix dimension and discard the size-0 ones.
   // Hopefully this can help with load balancing.
   bool isSorted = true;
   int maxSort = 10000;
   int count = 0;
   int temp;
   do {
      isSorted = true;
      for (int i=0; i < nTwoBodyChannels-1; i++) {
	 count++;
	 //
	 if (TwoBodyChannels[i].GetNumberKets() > TwoBodyChannels[i+1].GetNumberKets()) {
	    temp = SortedTwoBodyChannels[i];
	    SortedTwoBodyChannels[i] = SortedTwoBodyChannels[i+1];
	    SortedTwoBodyChannels[i+1] = temp;
	    isSorted = false;
	 }
      }
   } while (!isSorted and count < maxSort);
   //for (int i=nTwoBodyChannels-1; i >= 0; i--){
      //cout << "TwoBodyChannels[" << i << "].GetNumberKets()=" << TwoBodyChannels[i].GetNumberKets() << endl;
      //if (TwoBodyChannels[i].GetNumberKets() == 0) {
	// TwoBodyChannels.erase(TwoBodyChannels.begin() + i);
	 //continue;
      //}
   //}
	 
   //sort(
   //   SortedTwoBodyChannels.begin(),
   //   SortedTwoBodyChannels.end(),
   //   [this](int i, int j){
	// int in = TwoBodyChannels[i].GetNumberKets();
	// int jn = TwoBodyChannels[j].GetNumberKets();
        // return in > jn;
     // }
   //); // Neet to ensure GetNumberKets is handled properly?
   //cout << "Sorted TwoBodyChannels, moving to _CC." << endl;
   int temp2;
   count = 0;
   isSorted = true;
   do {
      isSorted = true;
      for (int i=0; i < nTwoBodyChannels-1; i++) {
	 count++;
	 //if (TwoBodyChannels_CC[i-1].GetNumberKets() == 0) {
	 //   TwoBodyChannels_CC.erase(TwoBodyChannels_CC.begin() + i-1);
	 //   continue;
	 //}
	 if (TwoBodyChannels_CC[i].GetNumberKets() > TwoBodyChannels_CC[i+1].GetNumberKets()) {
	    temp2 = SortedTwoBodyChannels_CC[i];
	    SortedTwoBodyChannels_CC[i] = SortedTwoBodyChannels_CC[i+1];
	    SortedTwoBodyChannels_CC[i+1] = temp2;
	    isSorted = false;
	 }
      }
   } while (!isSorted and count < maxSort);
   //for (int i=nTwoBodyChannels-1; i >= 0; i--){
      //cout << "TwoBodyChannels_CC[" << i << "].GetNumberKets()=" << TwoBodyChannels_CC[i].GetNumberKets() << endl;
      //if (TwoBodyChannels_CC[i].GetNumberKets() == 0) {
	// SortedTwoBodyChannels_CC.erase(TwoBodyChannels_CC.begin() + i);
	 //continue;
      //}
   //}
   //sort(
   //   SortedTwoBodyChannels_CC.begin(),
   //   SortedTwoBodyChannels_CC.end(),
   //   [this](int i, int j){ 
   //      int in = TwoBodyChannels[i].GetNumberKets();
//	 int jn = TwoBodyChannels[j].GetNumberKets();
   //      return in > jn;
   //   }  
   //); // Neet to ensure GetNumberKets is handled properly?
   //cout << "About to pop_back." << endl;
   while (  TwoBodyChannels[ SortedTwoBodyChannels.back() ].GetNumberKets() <1 ) SortedTwoBodyChannels.pop_back();
   while (  TwoBodyChannels_CC[ SortedTwoBodyChannels_CC.back() ].GetNumberKets() <1 ) SortedTwoBodyChannels_CC.pop_back();
   //for (int i=0; i < SortedTwoBodyChannels.size(); i++){
   //   cout << "Sorted TwoBodyChannels[" << i << "].GetNumberKets()=" << TwoBodyChannels[i].GetNumberKets() << endl;
   //}
   //for (int i=0; i < SortedTwoBodyChannels_CC.size(); i++){
   //   cout << "Sorted TwoBodyChannels_CC[" << i << "].GetNumberKets()=" << TwoBodyChannels_CC[i].GetNumberKets() << endl;
   //}
   //nTwoBodyChannels = TwoBodyChannels.size();
}


void ModelSpace::ClearVectors()
{
   holes.clear();         
   particles.clear();     
   core.clear();          
   valence.clear();       
   qspace.clear();        
   proton_orbits.clear();  
   neutron_orbits.clear();
   
   KetIndex_pp.clear();
   KetIndex_ph.clear();
   KetIndex_hh.clear();
   KetIndex_cc.clear();
   KetIndex_vc.clear();
   KetIndex_qc.clear();
   KetIndex_vv.clear();
   KetIndex_qv.clear();
   KetIndex_qq.clear();
   Ket_occ_hh.clear();
   Ket_occ_ph.clear();
   Ket_unocc_hh.clear();
   Ket_unocc_ph.clear();
   for (index_t Tz=0; Tz<3; ++Tz)
   {
     for (index_t parity=0;parity<2; ++parity)
     {
        MonopoleKets[Tz][parity].clear();
     }
   }

   Orbits.clear();
   Kets.clear();
   OneBodyChannels.clear();
   TwoBodyChannels.clear();
   TwoBodyChannels_CC.clear();
   SortedTwoBodyChannels.clear();
   SortedTwoBodyChannels_CC.clear();
}

void ModelSpace::SetSystemType(string str)
{
   SystemType = str;
}


double ModelSpace::GetSixJ(double j1, double j2, double j3, double J1, double J2, double J3)
{
// { j1 j2 j3 }
// { J1 J2 J3 }

//   unsigned long long int key = 20000000000*j1 + 200000000*j2 + 2000000*j3 + 20000*J1 + 200*J2 + 2*J3;
   unsigned long int key = (((unsigned long int) (2*j1)) << 30) +
                           (((unsigned long int) (2*j2)) << 24) +
                           (((unsigned long int) (2*j3)) << 18) +
                           (((unsigned long int) (2*J1)) << 12) +
                           (((unsigned long int) (2*J2)) <<  6) +
                            ((unsigned long int) (2*J3));

   auto it = SixJList.find(key);
   if (it != SixJList.end() ) return it->second;
   double sixj = AngMom::SixJ(j1,j2,j3,J1,J2,J3);
   #pragma omp critical
   SixJList[key] = sixj;
   return sixj;
}

void ModelSpace::PreCalculateMoshinsky() { PreCalculateMoshinsky( "harmonic" );}

void ModelSpace::PreCalculateMoshinsky( string basis )
{
//  if ( not MoshList.empty() ) return; // Already done calculated it...
  int Nmax = 0;
  cout << "systemBasis=" << basis << endl;
  if ( basis == "hydrogen" ) Nmax = 26; // just for testing;
  if ( basis == "harmonic" ) Nmax = E2max;
  #pragma omp parallel for schedule(dynamic,1)
  for (int N=0; N<=Nmax/2; ++N)
  {
   unordered_map<unsigned long long int,double> local_MoshList;
   for (int n=0; n<=min(N,Nmax/2-N); ++n)
   {
    for (int Lam=0; Lam<=Nmax-2*N-2*n; ++Lam)
    {
     int lam_max = (N==n ? min(Lam,Nmax-2*N-2*n-Lam) : Nmax-2*N-2*n-Lam);
     for (int lam=0; lam<=lam_max; ++lam)
     {
      int e2 = 2*N+Lam + 2*n+lam;
      for (int L=abs(Lam-lam); L<=Lam+lam; ++L)
      {
       for (int n1=0; n1<=N; ++n1)
       {
        for (int n2=0; n2<=min(n1,e2/2-n1); ++n2)
        {
         int l1max = n1==N? min(Lam,e2-2*n1-2*n2) : e2-2*n1-2*n2;
         for (int l1=0; l1<=l1max; ++l1 )
         {
          int l2 = e2-2*n1-2*n2-l1;
          if ( (l1+l2+lam+Lam)%2 >0 ) continue;
          if ( l2<abs(L-l1) or l2>L+l1 ) continue;
          // emax = 16, lmax = 32 -> good up to emax=32, which I'm nowhere near.
          unsigned long long int key =   ((unsigned long long int) N   << 40)
                                       + ((unsigned long long int) Lam << 34)
                                       + ((unsigned long long int) n   << 30)
                                       + ((unsigned long long int) lam << 26)
                                       + ((unsigned long long int) n1  << 22)
                                       + ((unsigned long long int) l1  << 16)
                                       + ((unsigned long long int) n2  << 12)
                                       + ((unsigned long long int) l2  << 6 )
                                       +  L;
          double mosh = AngMom::Moshinsky(N,Lam,n,lam,n1,l1,n2,l2,L);
          local_MoshList[key] = mosh;
         } // l1
        } // n2
       } // n1
      } // L
     } // lam
    } // Lam
   } // n
   #pragma omp critical
   MoshList.insert( local_MoshList.begin(), local_MoshList.end() );
  }
}

void ModelSpace::PreCalculateMoshinsky_FromList( vector<unsigned long long int>& mosh_list )
{
//  if ( not MoshList.empty() ) return; // Already done calculated it...
    bool gotZero=false;
    //#pragma omp parallel for
    for ( unsigned long long int it=0; it<mosh_list.size(); it++ )
    {
	//cout << "Calculating mosh for " << mosh_list[it] << endl;
	unsigned long long int temp;
	int N, Lam, n, lam, n1, l1, n2, l2, L;
	temp = mosh_list[it];
	if ( temp == 0 and gotZero == true ) continue;
	if ( temp == 0 and gotZero == false ) gotZero = true;
	L = temp%100;
	temp -= L;
	temp /=100;
	l2 = temp%100;
	temp -= l2;
	temp /=100;
	n2 = temp%100;
	temp -= n2;
	temp /=100;
	l1 = temp%100;
	temp -= l1;
	temp /=100;
	n1 = temp%100;
	temp -= n1;
	temp /=100;
	lam = temp%100;
	temp -= lam;
	temp /=100;
	n = temp%100;
	temp -= n;
	temp /=100;
	Lam = temp%100;
	temp -= Lam;
	temp /=100;
	N = temp%100;
	temp -= N;
	int phase_mosh = 1;
	int switches = 10;
	//cout << "N=" << N << " Lam=" << Lam << " n=" << n << " lam=" << lam << " n1=" << n1 << " l1=" << l1 << " n2=" << n2 << " L=" << L << endl;
	while (switches > 0)
	{
	    switches = 0;
	    if (n2>n1 or (n2==n1 and l2>l1))
   	    {
		swap(n1,n2);
		swap(l1,l2);
		phase_mosh *= phase(Lam+L);
		++switches;
	    }
	    if (n>N or (n==N and lam>Lam))
	    {
		swap(n,N);
		swap(lam,Lam);
		phase_mosh *= phase(l1 +L);
		++switches;
	    }
	    if (n1>N or (n1==N and l1>Lam) or (n1==N and l1==Lam and n2>n) or (n1==N and l1==Lam and n2==n and l2>lam) )
	    {
		swap(n1,N);
		swap(l1,Lam);
		swap(n2,n);
		swap(l2,lam);
		++switches;
//      phase_mosh *= phase(l2+lam); // This phase is given in Moshinsky and Brody, but with the current algorithm, it appears not to be required.
	    }
	}
          unsigned long long int key =   ((unsigned long long int) N   << 40)
                                       + ((unsigned long long int) Lam << 34)
                                       + ((unsigned long long int) n   << 30)
                                       + ((unsigned long long int) lam << 26)
                                       + ((unsigned long long int) n1  << 22)
                                       + ((unsigned long long int) l1  << 16)
                                       + ((unsigned long long int) n2  << 12)
                                       + ((unsigned long long int) l2  << 6 )
                                       +  L;
	auto iter = MoshList.find(key);
	/* unsigned long long int tkey = 0;
		tkey += pow(100,8)*N;
		tkey += pow(100,7)*Lam;
		tkey += pow(100,6)*n;
		tkey += pow(100,5)*lam;
		tkey += pow(100,4)*n1;
		tkey += pow(100,3)*l1;
		tkey += pow(100,2)*n2;
		tkey += 100*l2;
		tkey += L; */
   	if ( iter == MoshList.end() )
	{
	    //cout << "Making new moshinsky with tkey =" << tkey << endl;
            double mosh = AngMom::Moshinsky(N,Lam,n,lam,n1,l1,n2,l2,L);
	    #pragma omp critical
            MoshList[ key ] = mosh;
	}
    }
} 

double ModelSpace::GetMoshinsky( int N, int Lam, int n, int lam, int n1, int l1, int n2, int l2, int L)
{
  int phase_mosh = 1;
  int switches = 10;

  while (switches > 0)
  {
   switches = 0;
   if (n2>n1 or (n2==n1 and l2>l1))
   {
      swap(n1,n2);
      swap(l1,l2);
      phase_mosh *= phase(Lam+L);
      ++switches;
   }
   if (n>N or (n==N and lam>Lam))
   {
      swap(n,N);
      swap(lam,Lam);
      phase_mosh *= phase(l1 +L);
      ++switches;
   }

   if (n1>N or (n1==N and l1>Lam) or (n1==N and l1==Lam and n2>n) or (n1==N and l1==Lam and n2==n and l2>lam) )
   {
      swap(n1,N);
      swap(l1,Lam);
      swap(n2,n);
      swap(l2,lam);
      ++switches;
//      phase_mosh *= phase(l2+lam); // This phase is given in Moshinsky and Brody, but with the current algorithm, it appears not to be required.
   }
  }

          unsigned long long int key =   ((unsigned long long int) N   << 40)
                                       + ((unsigned long long int) Lam << 34)
                                       + ((unsigned long long int) n   << 30)
                                       + ((unsigned long long int) lam << 26)
                                       + ((unsigned long long int) n1  << 22)
                                       + ((unsigned long long int) l1  << 16)
                                       + ((unsigned long long int) n2  << 12)
                                       + ((unsigned long long int) l2  << 6 )
                                       +  L;

//   unsigned long long int key =  1000000000000 * N
//                                + 100000000000 * Lam
//                                +   1000000000 * n
//                                +    100000000 * lam
//                                +      1000000 * n1
//                                +       100000 * l1
//                                +         1000 * n2
//                                +          100 * l2
//                                +                 L;
   auto it = MoshList.find(key);
   if ( it != MoshList.end() )  return it->second * phase_mosh;
	/* unsigned long long int tkey = 0;
			tkey += pow(100,8)*N;
			tkey += pow(100,7)*Lam;
			tkey += pow(100,6)*n;
			tkey += pow(100,5)*lam;
			tkey += pow(100,4)*n1;
			tkey += pow(100,3)*l1;
	 		tkey += pow(100,2)*n2;
			tkey += 100*l2;
			tkey += L; */
   //cout << "Didn't find Moshinsky key, making a new one; tkey=" << tkey << endl;
   //#pragma omp critical

   //cout << "N=" << N << " Lam=" << Lam << " n=" << n << " lam=" << lam << " n1=" << n1 << " l1=" << l1 << " n2=" << n2 << " n2=" << n2 << " l2=" << l2 << endl;
   // if we didn't find it, we need to calculate it.
   double mosh = AngMom::Moshinsky(N,Lam,n,lam,n1,l1,n2,l2,L);
//   cout << "Shouldn't be here..." << N << " " << Lam << " " <<  n << " " << lam << " " << n1 << " " << l1 << " " << n2 << " " << l2 << " " << L << endl;
   //#pragma omp critial
   //{
       MoshList[key] = mosh;
   //}
   return mosh * phase_mosh;

}

struct my_f_params {int n; double l; int np; int Z;};

double
OsToHydroCoeff( double x, void * p )
{
	struct my_f_params * params = (struct my_f_params *)p;
        int n = (params->n);
        double l = (params->l);
	int np = (params->np);
	int Z = (params->Z);//2;//GetTargetZ(); //Fix later
	double c = Z / (n * BOHR_RADIUS); // 
	double m = 1; 			// Electron mass, in atomic units
	double h = 1; 			// Reduced plancks' constant, in atomic units
	double w = 13.605 * 2 / h; 	// wavelength of oscillator 13.605 *2 ? 
	double v = m*w/(2*h);
	       
	return pow( x, 2*l+2 ) * exp( -v*x*x - c*x ) * gsl_sf_laguerre_n(np, l+0.5, 2*v*x*x) * gsl_sf_laguerre_n(n-l-1, 2*l+1, 2*x*c);
}

void ModelSpace::GenerateOsToHydroCoeff(int nmax) {
    cout << "Entering GenerateOsToHydroCoeff." << endl;
    //OsToHydroCoeffList.resize(1000*46 + 10*nmax + 1*Lmax);

    int size = 1000;

    gsl_integration_workspace *work_ptr =
    gsl_integration_workspace_alloc (size);

    //gsl_set_error_handler_off();

    double lower_limit = 0;	/* start integral from lower_limit (to infinity) */
    double abs_error = 1.0e-6;	/* to avoid round-off problems */
    double rel_error = 1.0e-6;	/* the result will usually be much better */
    double result;		/* the result from the integration */
    double error;		/* the estimated error from the integration */

    gsl_function My_function;
    struct my_f_params alpha = {1,0,1,1};
    My_function.function = &OsToHydroCoeff;
    My_function.params = &alpha;
    int Z = GetTargetZ();
    int npmax = min(32, 4*Emax);

    //#pragma omp parallel for schedule(dynamic,1)
    for (int n = 1.; n <= nmax; n++)
    {
	for (int l = 0.; l < n and l <= Lmax; l++)
	{
	    double hydrogenCoeff = sqrt( pow(2*Z/(n * BOHR_RADIUS),3) * GetFactorial(n-l-1)/((2*n*GetFactorial(n+l)) ) ) * pow(2*Z/(n * BOHR_RADIUS),l);
	    //double temp = 0;
	    //#pragma parallel for schedule(dynamic,1)
	    for (int np = 0; np <= npmax; np++) // Should goto inf; throws NaN at np > 46; seems to throw at higher if you reduce errors
	    {
		alpha.n = n;
		alpha.l = l;
		alpha.np = np;
		alpha.Z = Z;
        	double OscilCoeff = sqrt( sqrt( 2* pow(13.605,3) / 3.14159 ) * pow(2, np+2*l+3) * GetFactorial(np) * pow(13.605,l) / gsl_sf_doublefact(2*np+2*l+1) ); // can split+cache
		
		//cout << "About to integrate; n=" << n << " l=" << l << " np=" << np << " hydrogenCoeff=" << hydrogenCoeff << " OscilCoeff=" << OscilCoeff << endl;

		gsl_integration_qagiu (&My_function,
					lower_limit,
					abs_error,
					rel_error,
					size,
					work_ptr,
					&result,
					&error);
	        //if ( isnan(result) == 1 ) continue;
		
		if ( std::isnan( result ) ) {
		    cout << "Hit NaN in Generating os Coeff." << endl;
		    continue;
		}
		int index = 1000*np + 10*n + l;
		cout << "Index =" << index << endl;
		cout << "Result=" << result << endl;
		//cout << "Hcoeff=" << hydrogenCoeff << endl;
		//cout << "Ocoeff=" << OscilCoeff << endl;
		//cout << "Coeff =" << OscilCoeff * hydrogenCoeff * result << endl;
		//if( hydrogenCoeff == 0 or result == 0 or OscilCoeff == 0 ){ // coeffs should never be zero, but just in case of rounding
		    //cout << "hydrogenCoeff=" << hydrogenCoeff << " temp=" << temp << endl;
		  //  OsToHydroCoeffList[index] = 0;
		//} else {
		    //cout << "Adding " << 1/(hydrogenCoeff*temp) << " to OsToHydro." << endl;
		if ( OsToHydroCoeffList[index] != 0 ) continue;
		OsToHydroCoeffList[index] = ( OscilCoeff * hydrogenCoeff * result ); // indexing should be good up to n = 9
		//cout << "OsToHydroCoeffList[index]=" << OsToHydroCoeffList[index] << endl;
	    	//}
	    	//temp += OscilCoeff*result;
	    } // np
	    //cout << "About to add for n=" << n << " l=" << l << endl;
	    
	    
	} // l
    } // n

    cout << "Exiting GenerateOsToHydroCoeff." << endl;
}

void ModelSpace::GenerateOsToHydroCoeff_fromlist( vector<int>& hy_list ) {
    cout << "Entering GenerateOsToHydroCoeff_fromList." << endl;
    //OsToHydroCoeffList.resize(1000*46 + 10*Emax + 1*Lmax); //Arbitrarily large size.

    int size = 1000;

    gsl_integration_workspace *work_ptr =
    gsl_integration_workspace_alloc (size);

    //gsl_set_error_handler_off();

    double lower_limit = 0;	/* start integral from lower_limit (to infinity) */
    double abs_error = 1.0e-6;	/* to avoid round-off problems */
    double rel_error = 1.0e-6;	/* the result will usually be much better */
    double result;		/* the result from the integration */
    double error;		/* the estimated error from the integration */

    gsl_function My_function;
    struct my_f_params alpha = {1,0,1,1};
    My_function.function = &OsToHydroCoeff;
    My_function.params = &alpha;
    int Z = GetTargetZ();

    //#pragma omp parallel for
    for ( unsigned long long int it=0; it<hy_list.size(); it++ )
    {
	auto iter = OsToHydroCoeffList.find(hy_list[it]);
	if (iter != OsToHydroCoeffList.end() ) continue;
	int temp = hy_list[it];
	if ( temp == 0 ) continue;
	alpha.l = temp%10;
	int l = alpha.l;
	temp -= alpha.l;
	temp /= 10;
	alpha.n = temp%100;
	int n = alpha.n;
	temp -= alpha.n;
	temp /= 100;
	alpha.np = temp;
	int np = alpha.np;
	alpha.Z = Z;
	double hydrogenCoeff = sqrt( pow(2*Z/(n * BOHR_RADIUS),3) * GetFactorial(n-l-1)/((2*n*GetFactorial(n+l)) ) ) * pow(2*Z/(n * BOHR_RADIUS),l);
        double OscilCoeff = sqrt( sqrt( 2* pow(13.605,3) / 3.14159 ) * pow(2, np+2*l+3) * GetFactorial(np) * pow(13.605,l) / gsl_sf_doublefact(2*np+2*l+1) ); // can split+cache
		
	//cout << "About to integrate; n=" << n << " l=" << l << " np=" << np << " hydrogenCoeff=" << hydrogenCoeff << " OscilCoeff=" << OscilCoeff << endl;
	gsl_integration_qagiu (&My_function,
				lower_limit,
				abs_error,
				rel_error,
				size,
				work_ptr,
				&result,
				&error);
        //if ( isnan(result) == 1 ) continue;
	
	if ( std::isnan( result ) ) {
	    cout << "Hit NaN in Generating os Coeff." << endl;
	    continue;
	}
	//int index = 1000*np + 10*n + l;
	//cout << "Index =" << index << endl;
	//cout << "Result=" << result << endl;


	#pragma omp critical
	OsToHydroCoeffList[hy_list[it]] = ( OscilCoeff * hydrogenCoeff * result ); // indexing should be good up to n = 99
    } // np
    //cout << "About to add for n=" << n << " l=" << l << endl;
    cout << "Exiting GenerateOsToHydroCoeff_fromList." << endl;
}

void ModelSpace::GenerateFactorialList(double m){
    cout << "Entering GenerateFactorialList for m=" << m << endl;
    factorialList.resize(m+1);
    //long double t=1;
    factorialList[0] = 1;
    #pragma omp parallel for
    for (int i=1; i < factorialList.size(); i++)
    {
	if (factorialList[i] == 0)
	    factorialList[i] = gsl_sf_fact(i);
	//t *= abs(i);
	//factorialList[i] = abs(t);
	//cout << "Generating factorial for " << i << "!=" << factorialList[i] << endl;
    }
}

double ModelSpace::GetFactorial(double m){
    #pragma omp critical
    if (m > factorialList.size())
    {
	cout << m << "! not in Factorial List, regenerating." << endl;
	GenerateFactorialList(m);
    }
    return factorialList[m];
}

void ModelSpace::PrecalculateNineJ( vector<unsigned long long int>& ninejList )
{
    //#pragma omp parallel
    for ( unsigned long long int it=0; it < ninejList.size(); it++ )
    {
	unsigned long long int temp;
	double j1, j2, J12, j3, j4, J34, J13, J24, J;
	temp = ninejList[it];
	J = temp%100;
	temp -= J;
	J /= 2;
	J24 = temp%100;
	temp -= J24;
	J24 /= 2;
	J13 = temp%100;
	temp -= J13;
	J13 /= 2;
	J34 = temp%100;
	temp -= J34;
	J34 /= 2;
	j4 = temp%100;
	temp -= j4;
	j4 /= 2;
	j3 = temp%100;
	temp -= j3;
	j3 /= 2;
	J12 = temp%100;
	temp -= J12;
	J12 /= 2;
	j2 = temp%100;
	temp -= j2;
	j2 /= 2;
	j1 = temp%100;
	j1 /= 2;

/*	int k1 = 2*j1;
	int k2 = 2*j2;
	int K12 = 2*J12;
	int k3 = 2*j3;
	int k4 = 2*j4;
	int K34 = 2*J34;
	int K13 = 2*J13;
	int K24 = 2*J24;
	int K = 2*J;
	array<int,9> klist = {k1,k2,K12,k3,k4,K34,K13,K24,K};
	array<double,9> jlist = {j1,j2,J12,j3,j4,J34,J13,J24,J};
	int imin = min_element(klist.begin(),klist.end()) - klist.begin();
	switch (imin)
	{
	   case 0:
		klist = {k4,K34,k3,K24,K,K13,k2,K12,k1};
		jlist = {j4,J34,j3,J24,J,J13,j2,J12,j1};
		break;
	   case 1:
		klist = {K13,K,K24,k3,K34,k4,k1,K12,k2};
		jlist = {J13,J,J24,j3,J34,j4,j1,J12,j2};
		break;
	   case 2:
		klist = {k3,k4,K34,K13,K24,K,k1,k2,K12};
		jlist = {j3,j4,J34,J13,J24,J,j1,j2,J12};
		break;
	   case 3:
		klist = {K12,k2,k1,K,K24,K13,K34,k4,k3};
		jlist = {J12,j2,j1,J,J24,J13,J34,j4,j3};
		break;
	   case 4:
		klist = {k1,K12,k2,K13,K,K24,k3,K34,k4};
		jlist = {j1,J12,j2,J13,J,J24,j3,J34,j4};
		break;
	   case 5:
		klist = {K13,K24,K,k1,k2,K12,k3,k4,K34};
		jlist = {J13,J24,J,j1,j2,J12,j3,j4,J34};
		break;
	   case 6:
		klist = {k2,K12,k1,k4,K34,k3,K24,K,K13};
		jlist = {j2,J12,j1,j4,J34,j3,J24,J,J13};
		break;
	   case 7:
		klist = {K12,k1,k2,K34,k3,k4,K,K13,K24};
		jlist = {J12,j1,j2,J34,j3,j4,J,J13,J24};
		break;
	   case 8:
		break;
   	}

	unsigned long long int key =   klist[0];
	unsigned long long int factor = 100;
	for (int i=1; i<9; ++i)
	{
	    key += klist[i]*factor;
	    factor *=100;
	} */
	//auto iter = NineJList.find(key);
	//if (iter == NineJList.end() )
	//{
	    //cout << "Calculating ninej with key=" << ninejList[it] << endl;
	    double ninej = AngMom::NineJ(j1,j2,J12,j3,j4,J34,J13,J24,J);
	    #pragma omp critical
	    NineJList[ninejList[it]] = ninej;
	//}
    }
}


double ModelSpace::GetNineJ(double j1, double j2, double J12, double j3, double j4, double J34, double J13, double J24, double J)
{
//   cout << "Calling GetNineJ" << endl;
   int k1 = 2*j1;
   int k2 = 2*j2;
   int K12 = 2*J12;
   int k3 = 2*j3;
   int k4 = 2*j4;
   int K34 = 2*J34;
   int K13 = 2*J13;
   int K24 = 2*J24;
   int K = 2*J;

   array<int,9> klist = {k1,k2,K12,k3,k4,K34,K13,K24,K};
   array<double,9> jlist = {j1,j2,J12,j3,j4,J34,J13,J24,J};
   int imin = min_element(klist.begin(),klist.end()) - klist.begin();
   switch (imin)
   {
      case 0:
       klist = {k4,K34,k3,K24,K,K13,k2,K12,k1};
       jlist = {j4,J34,j3,J24,J,J13,j2,J12,j1};
       break;
      case 1:
       klist = {K13,K,K24,k3,K34,k4,k1,K12,k2};
       jlist = {J13,J,J24,j3,J34,j4,j1,J12,j2};
       break;
      case 2:
       klist = {k3,k4,K34,K13,K24,K,k1,k2,K12};
       jlist = {j3,j4,J34,J13,J24,J,j1,j2,J12};
       break;
      case 3:
       klist = {K12,k2,k1,K,K24,K13,K34,k4,k3};
       jlist = {J12,j2,j1,J,J24,J13,J34,j4,j3};
       break;
      case 4:
       klist = {k1,K12,k2,K13,K,K24,k3,K34,k4};
       jlist = {j1,J12,j2,J13,J,J24,j3,J34,j4};
       break;
      case 5:
       klist = {K13,K24,K,k1,k2,K12,k3,k4,K34};
       jlist = {J13,J24,J,j1,j2,J12,j3,j4,J34};
       break;
      case 6:
       klist = {k2,K12,k1,k4,K34,k3,K24,K,K13};
       jlist = {j2,J12,j1,j4,J34,j3,J24,J,J13};
       break;
      case 7:
       klist = {K12,k1,k2,K34,k3,k4,K,K13,K24};
       jlist = {J12,j1,j2,J34,j3,j4,J,J13,J24};
       break;
      case 8:
       break;
   }

   unsigned long long int key =   klist[0];
   unsigned long long int factor = 100;
   for (int i=1; i<9; ++i)
   {
      key += klist[i]*factor;
      factor *=100;
   }
   auto it = NineJList.find(key);
   if (it != NineJList.end() )
   {
     return it->second;
   }
   //cout << "Missing NineJ, making a new one; key=" << key << endl;
   double ninej = AngMom::NineJ(jlist[0],jlist[1],jlist[2],jlist[3],jlist[4],jlist[5],jlist[6],jlist[7],jlist[8]);
   //cout << "NineJ calculated." << endl;
   //#pragma omp critical
   NineJList[key] = ninej;
   //cout << "Nine J added to list, returning." << endl;
   return ninej;

}


