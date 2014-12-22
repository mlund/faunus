#include <faunus/faunus.h>
using namespace Faunus;

typedef Space<Geometry::Cuboid,PointParticle> Tspace;

int main() {

  InputMap mcp( "gcmol.input" );

  Tspace spc( mcp );
  spc.load( "state", Tspace::RESIZE );

  auto pot = Energy::Nonbonded<Tspace, Potential::CoulombLJ>( mcp );
  Move::GreenGC<Tspace> gc( mcp, pot, spc );

  EnergyDrift sys;
  sys.init( Energy::systemEnergy( spc,pot,spc.p ) );

  cout << atom.info() << pot.info() << spc.info() << textio::header( "MC Simulation Begins!" );

  MCLoop loop( mcp );
  while ( loop[0] ) {
    while ( loop[1] ) {
      sys += gc.move();
    }
    cout << loop.timing();
  }

  sys.checkDrift( Energy::systemEnergy( spc,pot,spc.p ) );
  spc.save("state");

  UnitTest test( mcp );
  sys.test( test );
  gc.test( test );

  cout << loop.info() << sys.info() << spc.info() << gc.info() << test.info() << endl;

  return test.numFailed();
}

/**
 * @page example_GCMolecular Example: Grand Canonical Molecules
 *
 * This is a example of Grand Canonical Monte Carlo simulation
 * of rigid molecules.
 *
 * Run the example directly from the faunus directory:
 *
 *     $ make
 *     $ cd src/examples/
 *     $ ./gcmol.run
 *
 *
 * Input
 * =====
 * All listed files including the above C++ and python programs can be found in `src/examples/`
 *
 * gcmol.input
 * -----------
 * \include gcmol.input
 *
 * gcmol.json
 * ----------
 * State atom types
 * \include gcmol.json
 *
 * gcmol.cpp
 * ---------
 * \include gcmol.cpp
 */


