#include <faunus/faunus.h>
using namespace Faunus;

#ifdef CUBOID
typedef Geometry::Cuboid Tgeometry;   // specify geometry - here cube w. periodic boundaries
#else
typedef Geometry::Sphere Tgeometry;   // sphere with hard boundaries
#endif
typedef Potential::CoulombHS Tpairpot;// particle pair potential: primitive model

typedef Space<Tgeometry,PointParticle> Tspace;

int main() {
  cout << textio::splash();           // show faunus banner and credits
  
  InputMap mcp("polymers.json");      // open user input file
  MCLoop loop(mcp);                   // class for handling mc loops
  UnitTest test(mcp);                 // class for unit testing

  Tspace spc(mcp);                    // instantiate space and load molecules/atoms

  // Create Space and a Hamiltonian with three terms
  auto pot = Energy::Nonbonded<Tspace,Tpairpot>(mcp)
    + Energy::ExternalPressure<Tspace>(mcp) + Energy::Bonded<Tspace>();

  spc.load("state");                                // load old config. from disk (if any)

  Move::Propagator<Tspace> mv( mcp, pot, spc );
  Analysis::CombinedAnalysis analyze(mcp, pot, spc);

  cout << atom.info() << spc.info() << pot.info() << textio::header("MC Simulation Begins!");

  while ( loop[0] ) {  // Markov chain 
    while ( loop[1] ) {
      mv.move();
      analyze.sample();
    } // end of micro loop

    cout << loop.timing();

  } // end of macro loop

  mv.test(test);

  // print information
  cout << loop.info() << mv.info() << analyze.info() << test.info();

  return test.numFailed();
}
/** 
@page example_polymers Example: Polymers

This will simulate an arbitrary number of linear polymers in the NPT ensemble
with explicit salt particles and implicit solvent (dielectric continuum).
We include the following Monte Carlo moves:

- salt translation
- monomer translation
- polymer translation and rotation
- polymer crankshaft and pivot rotations
- isobaric volume move (NPT ensemble)

![Hardsphere polyelectrolytes with counter ions](polymers.png)

Run this example from the `examples` directory:

~~~~~~~~~~~~~~~~~~~
$ make
$ cd src/examples
$ python ./polymers.py
~~~~~~~~~~~~~~~~~~~

polymers.py
===========
This is a python run-script that generates all required
input as well as run the main program. This is also used
for running the unittest when running `make test`. 

@include polymers.py
 */

