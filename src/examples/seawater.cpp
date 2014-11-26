#include <faunus/faunus.h>
using namespace Faunus;                   // use Faunus namespace
typedef Space<Geometry::Cuboid> Tspace;   // Type of simulation space
typedef Potential::CoulombHS Tpair;       // and pair potential

int main() {
  ::atom.includefile("seawater.json");     // load atom properties
  InputMap in("seawater.input");           // open parameter file for user input
  Energy::Nonbonded<Tspace,Tpair> pot(in);// Hamiltonian, non-bonded only
  Tspace spc(in);                         // Simulation space, particles etc.
  Group salt;                             // Group for salt particles
  salt.addParticles(spc,in);              // Add according to user input
  Move::AtomicTranslation<Tspace> mv(in,pot,spc);// particle move class
  mv.setGroup(salt);                      // move class acts on salt group

  Analysis::WidomScaled<Tspace>
    widom( pot.pairpot.first.bjerrumLength(), 10 );
  widom.add( spc.p );

  MCLoop loop(in);
  while ( loop[0] )
    while ( loop[1] ) {
      mv.move( spc.p.size() );            // move salt randomly 100000 times
      widom.sample( spc.p, spc.geo );
    }

  cout << spc.info() + loop.info() + pot.info() + mv.info() + widom.info();
}
/**
 * @page example_seawater Example: Activity Coefficients in Seawater
 *
 * @includelineno seawater.cpp
 *
 * Run the code directly from the faunus directory:
 *
 *     $ make example_seawater
 *     $ cd src/examples/
 *     $ ./seawater.sh
 */
