#include <faunus/faunus.h>
#include <faunus/ewald.h>
using namespace Faunus;
using namespace Faunus::Potential;

//typedef CombinedPairPotential<HardSphereCap,CoulombGalore> Tpairpot; // pair potential
typedef HardSphereCap Tpairpot; // pair potential

typedef Geometry::Cuboid Tgeometry;   // geometry: cube w. periodic boundaries
typedef Space<Tgeometry,CapParticle> Tspace;

int main() {
  cout << textio::splash();           // show faunus banner and credits

  Tmjson mcp = openjson("capparticles.json"); // open JSON input file
  MCLoop loop(mcp);                   // class for handling mc loops

  Tspace spc(mcp);                    // simulation space

  auto pot = Energy::NonbondedVector<Tspace,Tpairpot>(mcp);

  spc.load("state");                  // load old config. from disk (if any)

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
  
  CapparticlePOVRay a;
  a.saveCapparticlePOVRay(mcp["analysis"]["povray"]["file"],spc);

  // perform unit tests (irrelevant for the simulation)
  UnitTest test(mcp);                    // class for unit testing
  mv.test(test);
  analyzer.test(test);

  // print information
  cout << loop.info() + mv.info() + analyzer.info() + test.info();

  return test.numFailed();
}

/**
  @page example_capparticles Example: Capped particles

  In this example we simulate 'Capparticles'. We only use a hard interaction.
  
  Information about the input file can be found in `src/examples/capparticles.py`.

  capparticles.cpp
  ========
  @includelineno examples/capparticles.cpp
*/
