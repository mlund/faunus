#include <faunus/faunus.h>
using namespace Faunus;                   // use Faunus namespace
typedef Space<Geometry::Cuboid> Tspace;   // Type of simulation space
typedef Potential::CoulombHS Tpair;       // and pair potential

int main() {
  InputMap in("seawater.json");           // open parameter file for user input
  Energy::Nonbonded<Tspace,Tpair> pot(in);// Hamiltonian, non-bonded only
  Tspace spc(in);                         // Simulation space, particles etc.

  Move::Propagator<Tspace> mv(in,pot,spc);// particle move class

  Analysis::CombinedAnalysis analyse(in,pot,spc);

  MCLoop loop(in);
  while ( loop[0] )
    while ( loop[1] ) {
      mv.move();
      analyse.sample();
    }

  cout << spc.info() + loop.info() + pot.info() + mv.info() + analyse.info();
}
/**
@page example_seawater Example: Activity Coefficients in Seawater

@includelineno seawater.cpp

Run the code directly from the faunus directory:

    $ make example_seawater
    $ cd src/examples/
    $ ./seawater.sh
 */
