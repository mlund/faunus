#include <faunus/faunus.h>
using namespace Faunus;
using namespace Faunus::Potential;

typedef Space<Geometry::Cuboid,PointParticle> Tspace;
typedef CombinedPairPotential<CoulombWolf,LennardJonesLB> Tpairpot;

int main() {

  InputMap mcp( "gcmol.input" );

  Tspace spc( mcp );
  spc.load( "state", Tspace::RESIZE );

  auto pot = Energy::Nonbonded<Tspace,Tpairpot>(mcp);
  Move::GreenGC<Tspace> gc( mcp, pot, spc );
  Move::TranslateRotate<Tspace> tr(mcp,pot,spc);

  EnergyDrift sys;
  sys.init( Energy::systemEnergy( spc,pot,spc.p ) );

  cout << atom.info() << pot.info() << textio::header( "MC Simulation Begins!" );

  MCLoop loop( mcp );
  while ( loop[0] ) {
    while ( loop[1] ) {
      if ( slump() > 0.5 ) 
        sys += gc.move();
      else 
        for ( auto cnt : spc.groupList() ) {
          auto it = slump.element( spc.groupList().begin(), spc.groupList().end() ) ;
          tr.setGroup( **it );
          sys += tr.move();
        }
    }
    cout << loop.timing();
  }

  sys.checkDrift( Energy::systemEnergy( spc,pot,spc.p ) );
  spc.save("state");
  FormatPQR::save("confout.pqr", spc.p, spc.geo.len);

  UnitTest test( mcp );
  sys.test( test );
  gc.test( test );
  tr.test( test );

  cout << loop.info() << sys.info() << spc.info() << tr.info() << gc.info() << test.info() << endl;

  return test.numFailed();
}

/**
 * @page example_GCMolecular Example: Grand Canonical Water
 *
 * This is a example of Grand Canonical Monte Carlo simulation
 * of rigid water molecules. The specified water activity corresponds
 * to the water vapour pressure in which with the liquid phase is in
 * equilibrium.
 *
 * ![Grand canonical water](NmuT.png)
 *
 * Run the example directly from the faunus directory:
 *
 *     $ make
 *     $ cd src/examples/
 *     $ ./gcmol.run
 *
 * Input
 * =====
 * All listed files including the above C++ and python programs can be found in `src/examples/`
 *
 * gcmol.input
 * -----------
 * @include gcmol.input
 *
 * gcmol.json
 * ----------
 * State atom types
 * @include gcmol.json
 *
 * gcmol.cpp
 * ---------
 * @include gcmol.cpp
 */


