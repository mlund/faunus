#include <faunus/faunus.h>
#include <faunus/ewald.h>
#define EWALD

using namespace Faunus;
using namespace Faunus::Potential;

typedef Space<Geometry::Cuboid> Tspace;
#ifdef EWALD
typedef LennardJonesLB Tpairpot;
#else
typedef CombinedPairPotential<CoulombGalore,LennardJonesLB> Tpairpot;
#endif

int main() {

  cout << textio::splash();      // show faunus banner and credits
  Tmjson mcp = openjson("water2.json");// read input file
  Tspace spc(mcp);
  
  // Energy functions and space
#ifdef EWALD
  auto pot = Energy::NonbondedEwald<Tspace,Tpairpot>(mcp)
    + Energy::ExternalPressure<Tspace>(mcp);
#else
  auto pot = Energy::NonbondedCutg2g<Tspace,Tpairpot>(mcp)
    + Energy::ExternalPressure<Tspace>(mcp);
#endif

  spc.load("state"); // load old config. from disk (if any)

  // Markov moves and analysis
  Analysis::CombinedAnalysis analyzer(mcp,pot,spc);
  Move::Propagator<Tspace> mv( mcp, pot, spc );

  cout << atom.info() + spc.info() + pot.info() + textio::header("MC Simulation Begins!");
  
  MCLoop loop(mcp);    // class for handling mc loops
  while ( loop[0] ) {          // Markov chain 
    while ( loop[1] ) {
      mv.move();
      analyzer.sample();
    } // end of micro loop

    cout << loop.timing();
  } // end of macro loop

  // perform unit 
  UnitTest test(mcp);
  mv.test(test);

  // print information
  cout << loop.info() + mv.info() + analyzer.info() + test.info();

  return test.numFailed();
}

/**  @page example_water2 Example: SPC/E Water (V2)

 This will simulate SPC/E water in a cubic volume using
 the q-potential for electrostatic interactions. When
 `order` in input is set to 1 then the q-potential is
 equal to the Wolf method. This version uses a fake 
 cell list to discard interactions beyond a specified 
 water-water mass-center cutoff.
 
 The simulation results can be directly compared with
 results presented in DOI:10.1063/1.476482.

 Run this example from the main faunus directory:

 ~~~~~~~~~~~~~~~~~~~
 $ make example_water2
 $ cd src/examples
 $ ./water2.py
 ~~~~~~~~~~~~~~~~~~~

 ![Water](water.png)

 water2.cpp
 ==========

 @includelineno examples/water2.cpp

 water2.json
 -----------

 @includelineno examples/water2.json

*/
