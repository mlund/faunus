#include <faunus/faunus.h>
#include <faunus/multipole.h>
#include <faunus/ewald.h>

//#define EWALD

using namespace Faunus;
using namespace Faunus::Move;
using namespace Faunus::Potential;

typedef Space<Geometry::Cuboid,DipoleParticle> Tspace; 
#ifdef EWALD
typedef LennardJonesLB Tpair;
#else
typedef CombinedPairPotential<LennardJonesLB,DipoleDipole> Tpair;
#endif

int main() {
  InputMap in("stockmayer.json");                // open parameter file for user input
  Tspace spc(in);                                // sim.space, particles etc.

#ifdef EWALD
    Energy::NonbondedEwald<Tspace,Tpair> pot(in); // non-bonded only
#else
    Energy::NonbondedVector<Tspace,Tpair> pot(in); // non-bonded only
#endif

#ifdef __POLARIZE
  Move::Propagator<Tspace,true> mv(in,pot,spc);
#else
  Move::Propagator<Tspace,false> mv(in,pot,spc);
#endif

  spc.load("state");

  Analysis::DipoleAnalysis dian(spc,in);
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
/**  @page example_stockmayer Example: Stockmayer potential

  This will simulate a Stockmayer potential in a cubic box.

  Run this example from the `examples` directory:

  ~~~~~~~~~~~~~~~~~~~
  $ make
  $ cd src/examples
  $ ./stockmayer.run
  ~~~~~~~~~~~~~~~~~~~

  stockmayer.cpp
  ==============

  @includelineno examples/stockmayer.cpp

*/
