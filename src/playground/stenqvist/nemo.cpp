#include <faunus/faunus.h>
#include <faunus/ewald.h>
#include "faunus/multipole.h"

using namespace Faunus;
using namespace Faunus::Potential;

typedef Space<Geometry::Cuboid, DipoleParticle> Tspace;
#ifdef EWALD
typedef LennardJonesLB Tpairpot;
#else
typedef IonIonQ Tpairpot0;
typedef CombinedPairPotential<LennardJonesLB,Tpairpot0> Tpairpot;
#endif

template<class Tpairpot, class Tid>
bool savePotential(Tpairpot pot, Tid ida, Tid idb, string file) {
  std::ofstream f(file.c_str());
  if (f) {
    DipoleParticle a,b;
    a=atom[ida];
    b=atom[idb];
    a.mu = Point(1,0,0);
    b.mu = Point(1,0,0);
    for (double r=0.5; r<=4.5; r+=0.05) {
      f << std::left << std::setw(10) << r << " "
        << pot(a,b,Point(r,0,0)) << endl;
    }
    return true;
  }
  return false;
}

int main() {
  InputMap mcp("nemo.json");  

#ifdef EWALD
  auto pot = Energy::NonbondedEwald<Tspace,Tpairpot>(mcp)
    + Energy::ExternalPressure<Tspace>(mcp);
#else
  auto pot = Energy::NonbondedVector<Tspace,Tpairpot>(mcp)
    + Energy::ExternalPressure<Tspace>(mcp);
#endif
  auto pot0 = Energy::NonbondedVector<Tspace,LennardJonesLB>(mcp);
  auto pot1 = Energy::ExternalPressure<Tspace>(mcp);

  Tspace spc(mcp);
  spc.load("state");
  auto waters = spc.findMolecules("water");
  
  vector<double> energy_T;
  vector<double> energy_LJ;
  vector<double> energy_EP;

  Move::Propagator<Tspace> mv( mcp, pot, spc );
  Analysis::RadialDistribution<> rdfOO(0.05);
  Analysis::RadialDistribution<> rdfOH(0.05);
  Analysis::RadialDistribution<> rdfHH(0.05);
  
  savePotential(Tpairpot(mcp), atom["OW"].id, atom["OW"].id, "potential.dat");

  EnergyDrift sys; 
  sys.init( Energy::systemEnergy(spc,pot,spc.p)  ); 
  Analysis::MultipoleAnalysis dian(spc,mcp);

  MCLoop loop(mcp);   
  while ( loop[0] ) {       
    while ( loop[1] ) {

      sys+=mv.move();
      dian.samplePP(spc);

      double rnd = slump();
      if ( rnd>0.9 ) {
        rdfOO.sample( spc, atom["OW"].id, atom["OW"].id ); // O-O g(r)
        rdfOH.sample( spc, atom["OW"].id, atom["HW"].id ); // O-H g(r)
        rdfHH.sample( spc, atom["HW"].id, atom["HW"].id ); // H-H g(r)
        
        energy_T.push_back(Energy::systemEnergy(spc,pot,spc.p));
	energy_LJ.push_back(Energy::systemEnergy(spc,pot,spc.p));
	energy_EP.push_back(Energy::systemEnergy(spc,pot,spc.p));
      }
    } // end of micro loop

    sys.checkDrift(Energy::systemEnergy(spc,pot,spc.p)); 
    cout << loop.timing();
  } // end of macro loop

  spc.save("state");
  rdfOO.save("rdfOO.dat");
  rdfOH.save("rdfOH.dat");
  rdfHH.save("rdfHH.dat");
  dian.save();

  cout << spc.info() + pot.info() + mv.info() + sys.info() + dian.info();
  
  
  string file = "energy.dat";
  std::ofstream f(file.c_str());
  if (f) {
    for (int i = 0; i < energy_T.size(); i++)
      f << std::left << std::setw(10) << energy_T.at(i) << "     " << energy_LJ.at(i) << "      " << energy_EP.at(i) << endl;
    return true;
  }

  return 0;
}