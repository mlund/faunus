#include <faunus/faunus.h>
#include <faunus/ewald.h>

using namespace Faunus;
using namespace Faunus::Potential;

//#define EWALD

typedef Space<Geometry::Cuboid> Tspace;
#ifdef EWALD
typedef LennardJonesLB_SF Tpairpot;
#else
typedef CombinedPairPotential<CoulombWolf,LennardJonesLB> Tpairpot;
#endif

int main() {

  cout << textio::splash();      // show faunus banner and credits
  InputMap mcp("water2.json");   // read input file
  Tspace spc(mcp);
  spc.load("state"); // load old config. from disk (if any)

  // Energy functions and space
#ifdef EWALD
  auto pot = Energy::NonbondedEwald<Tspace,Tpairpot>(mcp)
    + Energy::ExternalPressure<Tspace>(mcp);
#else
  auto pot = Energy::NonbondedCutg2g<Tspace,Tpairpot>(mcp)
    + Energy::ExternalPressure<Tspace>(mcp);
#endif
  
  auto waters = spc.findMolecules("water");

  // Markov moves and analysis
  Move::Propagator<Tspace> mv( mcp, pot, spc );
  Analysis::RadialDistribution<> rdf(0.05);
  Analysis::MultipoleAnalysis<Tspace> dian(spc,pot,mcp);
  FormatXTC xtc(1000);

  EnergyDrift sys;   // class for tracking system energy drifts
#ifdef EWALD
  pot.setSpace(spc);
#endif

  sys.init( Energy::systemEnergy(spc,pot,spc.p)  ); // store total energy
  cout << atom.info() + spc.info() + pot.info() + textio::header("MC Simulation Begins!");
  
  MCLoop loop(mcp);    // class for handling mc loops
  while ( loop[0] ) {          // Markov chain 
    while ( loop[1] ) {

      sys+=mv.move();

      double rnd = slump();
      if ( rnd>0.9 ) {
        rdf.sample( spc, atom["OW"].id, atom["OW"].id ); // O-O g(r)
        if ( rnd > 0.99 ) {
          xtc.setbox( spc.geo.len );
          xtc.save( "traj.xtc", spc.p );
        }
        dian.sample(spc);
      }
    } // end of micro loop

    sys.checkDrift(Energy::systemEnergy(spc,pot,spc.p)); // energy drift?

    cout << loop.timing();
  } // end of macro loop

  // save to disk
  FormatPQR::save("confout.pqr", spc.p, spc.geo.len);
  spc.save("state");
  rdf.save("rdf.dat");

  // perform unit 
  UnitTest test(mcp);
  sys.test(test);
  mv.test(test);

  // print information
  cout << loop.info() + sys.info() + mv.info() + pot.info() + dian.info() + test.info();

  return test.numFailed();
}

/**  @page example_water2 Example: SPC Water (V2)

 This will simulate SPC water in a cubic volume using
 the Wolf method for electrostatic interactions.
 This version uses a fake cell list to discard
 interactions beyond a specified water-water mass-center
 cutoff.

 Run this example from the main faunus directory:

 ~~~~~~~~~~~~~~~~~~~
 $ make example_water2
 $ cd src/examples
 $ ./water2.run
 ~~~~~~~~~~~~~~~~~~~

 ![Water](water.png)

 water2.cpp
 ==========

 @includelineno examples/water2.cpp

 water2.json
 -----------

 @includelineno examples/water2.json

*/
