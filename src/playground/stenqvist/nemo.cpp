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
  ::atom.includefile("stockmayer.json");         // load atom properties
  InputMap in("stockmayer.input");               // open parameter file for user input
  Energy::NonbondedVector<Tspace,Tpair> pot(in); // non-bonded only
  EnergyDrift sys;                               // class for tracking system energy drifts
  Tspace spc(in);                // create simulation space, particles etc.
  Group sol;
  sol.addParticles(spc, in);                     // group for particles
  MCLoop loop(in);                               // class for mc loop counting
  Analysis::RadialDistribution<> rdf(0.05);       // particle-particle g(r)
  Analysis::Table2D<double,Average<double> > mucorr(0.1);       // particle-particle g(r)
  TmoveTran trans(in,pot,spc);
  TmoveRot rot(in,pot,spc);
  trans.setGroup(sol);                                // tells move class to act on sol group
  rot.setGroup(sol);                                  // tells move class to act on sol group
  spc.load("state");
  spc.p = spc.trial;
  UnitTest test(in);               // class for unit testing
  Analysis::DielectricConstant gdc(spc);
  FormatXTC xtc(spc.geo.len.norm());
  
  sys.init( Energy::systemEnergy(spc,pot,spc.p)  );   // initial energy
  while ( loop.macroCnt() ) {                         // Markov chain 
    while ( loop.microCnt() ) {
      if (slp_global() > 0.5)
        sys+=trans.move( sol.size() );                // translate
      else
        sys+=rot.move( sol.size() );                  // rotate
        
      gdc.sampleDP(spc);
      gdc.info();
      if (slp_global()<1.5)
        for (auto i=sol.front(); i<sol.back(); i++) { // salt rdf
          for (auto j=i+1; j<=sol.back(); j++) {
            double r =spc.geo.dist(spc.p[i],spc.p[j]); 
            rdf(r)++;
            mucorr(r) += spc.p[i].mu.dot(spc.p[j].mu);
          }
        }
      if (slp_global()>0.99)
        xtc.save(textio::prefix+"out.xtc", spc.p);  
    }    
    cout << gdc.info() << endl;
    //gdc.getKirkwoodFactor();
    sys.checkDrift(Energy::systemEnergy(spc,pot,spc.p)); // compare energy sum with current
    cout << loop.timing();
  }
  
  // perform unit tests
  trans.test(test);
  rot.test(test);
  sys.test(test);

  FormatPQR().save("confout.pqr", spc.p);
  rdf.save("gofr.dat");                               // save g(r) to disk
  mucorr.save("mucorr.dat");                               // save g(r) to disk
  std::cout << spc.info() + pot.info() + trans.info()
    + rot.info() + sys.info() + test.info(); // final info
  spc.save("state");
  gdc.kusalik(spc);
  
  return test.numFailed();
}
/**  @page example_stockmayer Example: Stockmayer potential
 *
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