#include <faunus/faunus.h>
using namespace Faunus;                               // use Faunus namespace
typedef Geometry::Cuboid Tgeo;                        // select simulation geometry and pair potential
typedef Potential::CombinedPairPotential<Potential::Coulomb, Potential::LennardJones> Tpair;

int main() {
  atom.includefile("atomlist.in");                    // load atom properties
  InputMap in("parameters.in");                       // open parameter file for user input
  Energy::Hamiltonian pot;                            // Hamiltonian - defines the energy field
  pot.create( Energy::Nonbonded<Tpair,Tgeo>(in) );    // add energy term for non-bonded interactions
  Space spc( pot.getGeometry() );                     // create simulation space, particles etc.
  GroupAtomic salt(spc, in);                          // group for salt particles
  Move::AtomicTranslation mv(in, pot, spc);           // particle move class
  mv.setGroup(salt);                                  // tells move class to act on salt group
  mv.move(1000);                                      // move salt randomly 1000 times
  std::cout << pot.info() << mv.info();               // final information
}
