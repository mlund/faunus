#include <faunus/faunus.h>
#include <faunus/ewald.h>

using namespace Faunus;
using namespace Faunus::Potential;

typedef Space<Geometry::Cuboid> Tspace;
#ifdef EWALD
typedef LennardJonesLB Tpairpot;
#else
typedef CombinedPairPotential<CoulombWolf,LennardJonesLB> Tpairpot;
#endif

int main() {

  cout << textio::splash();      // show faunus banner and credits
  InputMap mcp("water2.input");  // read input file

  // Energy functions and space
#ifdef EWALD
  auto pot = Energy::NonbondedEwald<Tspace,Tpairpot>(mcp)
    + Energy::ExternalPressure<Tspace>(mcp);
#else
  auto pot = Energy::NonbondedCutg2g<Tspace,Tpairpot>(mcp)
    + Energy::ExternalPressure<Tspace>(mcp);
#endif

  Tspace spc(mcp);
  spc.load("state"); // load old config. from disk (if any)
  auto waters = spc.findMolecules("water");

  // Markov moves and analysis
  Move::TranslateRotate<Tspace> gmv(mcp,pot,spc);
  Move::Isobaric<Tspace> iso(mcp,pot,spc);
  Analysis::RadialDistribution<> rdf(0.05);

  EnergyDrift sys;   // class for tracking system energy drifts
  sys.init( Energy::systemEnergy(spc,pot,spc.p)  ); // store total energy

  cout << atom.info() + spc.info() + pot.info() + textio::header("MC Simulation Begins!");

  MCLoop loop(mcp);    // class for handling mc loops
  while ( loop[0] ) {          // Markov chain 
    while ( loop[1] ) {
      int i=slump.range(0,1);
      int k=waters.size();
      Group g;
      switch (i) {
        case 0:
          while (k-->0) {
            gmv.setGroup( **slump.element(waters.begin(), waters.end()) );
            sys+=gmv.move();   // translate/rotate polymers
          }
          break;
        case 1:
          sys+=iso.move();
          break;
      }

      // sample oxygen-oxygen rdf
      if (slump()>0.9) {
        auto id = atom["OW"].id;
        rdf.sample(spc,id,id);
      }

    } // end of micro loop

    sys.checkDrift(Energy::systemEnergy(spc,pot,spc.p)); // energy drift?

    cout << loop.timing();
  } // end of macro loop

  // save to disk
  FormatPQR::save("confout.pqr", spc.p);
  spc.save("state");
  rdf.save("rdf.dat");
  spc.save("state");
  FormatPQR::save("confout.pqr", spc.p, spc.geo.len);

  // perform unit 
  UnitTest test(mcp);
  iso.test(test);
  gmv.test(test);
  sys.test(test);

  // print information
  cout << loop.info() + sys.info() + gmv.info() + iso.info() + test.info();

  return test.numFailed();
}

/**  @page example_water2 Example: SPC Water (V2)
 *
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
 ============

 @includelineno examples/water2.cpp

*/
