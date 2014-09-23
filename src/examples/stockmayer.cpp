#include <faunus/faunus.h>
#include <faunus/multipole.h>
using namespace Faunus;                          // use Faunus namespace
using namespace Faunus::Move;
using namespace Faunus::Potential;

typedef Space<Geometry::Cuboid,DipoleParticle> Tspace; 
typedef CombinedPairPotential<LennardJones,DipoleDipole> Tpair;

#ifdef POLARIZE
typedef Move::PolarizeMove<AtomicTranslation<Tspace> > TmoveTran;
typedef Move::PolarizeMove<AtomicRotation<Tspace> > TmoveRot;
#else
typedef Move::AtomicTranslation<Tspace> TmoveTran;   
typedef Move::AtomicRotation<Tspace> TmoveRot;
#endif

int main() {
  //::atom.includefile("stockmayer.json");         // load atom properties
  InputMap in("stockmayer.input");               // open parameter file for user input
  Energy::NonbondedVector<Tspace,Tpair> pot(in); // non-bonded only
  EnergyDrift sys;                               // class for tracking system energy drifts
  Tspace spc(in);                // create simulation space, particles etc.
  Group sol;
  sol.addParticles(spc, in);                     // group for particles
  MCLoop loop(in);                               // class for mc loop counting
  TmoveTran trans(in,pot,spc);
  TmoveRot rot(in,pot,spc);
  trans.setGroup(sol);                                // tells move class to act on sol group
  rot.setGroup(sol);                                  // tells move class to act on sol group
  spc.load("state");
  spc.trial = spc.p;
  UnitTest test(in);               // class for unit testing
  Analysis::DipoleAnalysis dian(spc,in);
  DipoleWRL sdp;
  FormatXTC xtc(spc.geo.len.norm());
  sys.init( Energy::systemEnergy(spc,pot,spc.p)  );   // initial energy
  while ( loop[0] ) {                         // Markov chain 
    while ( loop[1] ) {
      if (slp_global() > 0.5)
        sys+=trans.move( sol.size() );                // translate
      else
        sys+=rot.move( sol.size() );                  // rotate

      if (slp_global()<0.5)
        dian.sampleMuCorrelationAndKirkwood(spc);
      if (slp_global()>0.99)
        xtc.save(textio::prefix+"out.xtc", spc.p);  
      dian.sampleDP(spc);
    }    
    sys.checkDrift(Energy::systemEnergy(spc,pot,spc.p)); // compare energy sum with current
    cout << loop.timing() << std::flush;
  }

  // perform unit tests
  trans.test(test);
  rot.test(test);
  sys.test(test);
  sdp.saveDipoleWRL("stockmayer.wrl",spc,sol);
  FormatPQR().save("confout.pqr", spc.p);
  std::cout << spc.info() + pot.info() + trans.info()
    + rot.info() + sys.info() + test.info() + dian.info(); // final info
  dian.save();
  spc.save("state");

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
 ============

 @includelineno examples/stockmayer.cpp

*/
