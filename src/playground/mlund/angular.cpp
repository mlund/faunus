#include<faunus/faunus.h>

using namespace Faunus;                   // use Faunus namespace
typedef Space<Geometry::Cuboid> Tspace;   // Type of simulation space
typedef Potential::CoulombLJ Tpair;       // and pair potential

int main() {
  ::atom.includefile("../../examples/minimal.json");     // load atom properties
  InputMap in("../../examples/minimal.input");           // open parameter file for user input
  Energy::Nonbonded<Tspace,Tpair> pot(in);// Hamiltonian, non-bonded only
  Tspace spc(in);                         // Simulation space, particles etc.
  Group salt;                             // Group for salt particles
  salt.addParticles(spc,in);              // Add according to user input
  Move::AtomicTranslation<Tspace> mv(in,pot,spc);// particle move class
  mv.setGroup(salt);                      // move class acts on salt group
  //mv.move(1e5);                           // move salt randomly 100000 times
  cout << spc.info() + pot.info() + mv.info(); // final information

  Energy::Manybody<Tspace> three(spc);
  three.add( Potential::Angular( {0,1,2}, 30., 0.5 ) );
  cout << three.external(spc.p) << endl;
  cout << three.info();
}
