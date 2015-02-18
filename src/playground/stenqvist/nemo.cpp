#define DIPOLEPARTICLE
#include <faunus/faunus.h>
#include <faunus/multipole.h>
#include <functional>
#include <iostream>
using namespace Faunus;                     
using namespace Faunus::Move;
using namespace Faunus::Potential;

typedef Space<Geometry::Cuboid,DipoleParticle> Tspace; 
typedef DipoleDipole TpairDD;
typedef LennardJonesLB TpairLJ;
typedef CombinedPairPotential<TpairLJ,TpairDD> Tpair;

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
  InputMap in("stockmayer.json");                // open parameter file for user input
  Tspace spc(in);                                // sim.space, particles etc.

  Energy::NonbondedVector<Tspace,Tpair> pot(in); // non-bonded only
  //savePotential(pot, "sol", "sol","potential.dat");

#ifdef __POLARIZE
  Move::Propagator<Tspace,true> mv(in,pot,spc);
#else
  Move::Propagator<Tspace,false> mv(in,pot,spc);
#endif

  spc.load("state");

  Analysis::MultipoleAnalysis dian(spc,in);
  DipoleWRL sdp;
  FormatXTC xtc(spc.geo.len.norm());

  EnergyDrift sys;                               // class for tracking system energy drifts
  sys.init( Energy::systemEnergy(spc,pot,spc.p)  );// initial energy

  MCLoop loop(in);                               // class for mc loop counting
  while ( loop[0] ) {                            // Markov chain 
    while ( loop[1] ) {
      sys += mv.move();
      if (slump()>0.5)
        dian.sampleMuCorrelationAndKirkwood(spc);
      if (slump()>0.99)
        xtc.save(textio::prefix+"out.xtc", spc.p);  
      dian.sampleDP(spc);
    }    
    sys.checkDrift(Energy::systemEnergy(spc,pot,spc.p)); // compare energy sum with current
    cout << loop.timing() << std::flush;
  }

  UnitTest test(in);
  mv.test(test);
  sys.test(test);
  sdp.saveDipoleWRL("stockmayer.wrl", spc, Group(0, spc.p.size()-1) );
  FormatPQR().save("confout.pqr", spc.p);
  dian.save();
  spc.save("state");

  std::cout << spc.info() + pot.info() + mv.info()
    + sys.info() + test.info() + dian.info(); // final info

  return test.numFailed();
}
