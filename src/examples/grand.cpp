#include <faunus/faunus.h>
#include <faunus/ewald.h>
using namespace Faunus;
using namespace Faunus::Potential;

typedef Geometry::Sphere Tgeometry;             // cube w. periodic boundaries
typedef Space<Tgeometry,PointParticle> Tspace;

int main() {
  InputMap mcp("grand.input");                  // open user input file

  Tspace spc(mcp);                              // simulation space
  Energy::Nonbonded<Tspace,CoulombHS> pot(mcp); // hamiltonian
  pot.setSpace(spc);                            // share space w. hamiltonian

  Group salt;                                   // group for salt particles
  salt.addParticles(spc, mcp);
  spc.load("state",Tspace::RESIZE);             // load old config. from disk (if any)

  // Two different Widom analysis methods
  double lB = Coulomb(mcp).bjerrumLength();     // get bjerrum length
  Analysis::Widom<PointParticle> widom1;        // widom analysis (I)
  Analysis::WidomScaled<PointParticle> widom2(lB,10);// ...and (II)
  widom1.add(spc.p);
  widom2.add(spc.p);

  Move::GrandCanonicalSalt<Tspace> gc(mcp,pot,spc,salt);
  Move::AtomicTranslation<Tspace> mv(mcp,pot,spc);
  mv.setGroup(salt);

  EnergyDrift sys;                              // class for tracking system energy drifts
  sys.init(Energy::systemEnergy(spc,pot,spc.p));// store initial total system energy

  cout << atom.info() + spc.info() + pot.info() + "\n";

  MCLoop loop(mcp);                             // class for handling mc loops
  while ( loop.macroCnt() ) {
    while ( loop.microCnt() ) {
      if (slp_global() < 0.5)
        sys+=mv.move( salt.size() );            // translate salt
      else 
        sys+=gc.move();                         // grand canonical exchange

      widom1.sample(spc,pot,10);
      widom2.sample(spc.p,spc.geo);
    }                                           // end of micro loop

    sys.checkDrift(Energy::systemEnergy(spc,pot,spc.p)); // calc. energy drift
    cout << loop.timing();
  }                                             // end of macro loop

  FormatPQR::save("confout.pqr", spc.p);        // PQR snapshot for VMD etc.
  spc.save("state");                            // final simulation state

  UnitTest test(mcp);                           // class for unit testing
  gc.test(test);
  mv.test(test);
  sys.test(test);

  cout << loop.info() + sys.info() + mv.info() + gc.info() + test.info()
    + widom1.info() + widom2.info();

  return test.numFailed();
}
/**
  @page example_grand Example: Grand Canonical Salt

  This is an example of a grand caninical salt solution (NmuT).

  We have the following MC moves:
  - salt translation
  - salt exchange with virtual bulk - i.e. constant chemical potential

  Information about the input file can be found in `src/examples/grand.run`.

  grand.cpp
  ========
  \includelineno examples/grand.cpp
*/
