#include <faunus/faunus.h>
using namespace Faunus;
using namespace Faunus::Potential;

typedef CombinedPairPotential<Coulomb,HardSphere> Tpairpot; // pair potential

typedef Geometry::SphereSurface Tgeometry;   // geometry: spherical-surface
typedef Space<Tgeometry,PointParticle> Tspace;

int main() {
  cout << textio::splash();           // show faunus banner and credits

  Tmjson mcp = openjson("spheresurface.json"); // open JSON input file
  MCLoop loop(mcp);                   // class for handling mc loops
  
  Tspace spc(mcp);

  auto pot = Energy::Nonbonded<Tspace,Tpairpot>(mcp);
  
  spc.load("state");                               // load old config. from disk (if any)
  
  Analysis::CombinedAnalysis analyzer(mcp,pot,spc);
  
  Move::Propagator<Tspace> mv(mcp,pot,spc);

  cout << atom.info() + spc.info() + pot.info() + textio::header("MC Simulation Begins!");

  while ( loop[0] ) {  // Markov chain 
    while ( loop[1] ) {
      mv.move();
      analyzer.sample();
    } // end of micro loop

    cout << loop.timing();

  } // end of macro loop

  // perform unit tests (irrelevant for the simulation)
  UnitTest test(mcp);                    // class for unit testing
  mv.test(test);
  analyzer.test(test);
  
  // print information
  cout << loop.info() + mv.info() + analyzer.info() + test.info();

  return test.numFailed();
}

/**
  @page example_spheresurface Example: Charges on a 2D hypersphere-surface

  The `spheresurface.cpp` program can be used to get evenly spread particles on a surface of a sphere.

  We have the following MC moves:
  - particle translation

  Information about the input file can be found in `src/examples/spheresurface.py`.

  spheresurface.json
  =========
  @includelineno examples/spheresurface.json

  spheresurface.cpp
  ========
  @includelineno examples/spheresurface.cpp

*/
