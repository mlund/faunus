/*!
 * \page example_minimal Example: Minimal
 * This is a minimal example of how to set up a simulation within the faunus
 * framework. Here we simulate charged Lennard-Jones particles in a cubic simulation
 * box with periodic boundaries and the minimum image convention.
 * Two input/data files are required and we here show the implementation in both C++
 * and Python.
 * \section minimal_cpp minimal.cpp
 * \include minimal.cpp
 * \section minimal_py minimal.py
 * \include minimal.py
 * \section minimal_input Input
 * \subsection minimal_atoms minimal.atoms
 * \include minimal.atoms
 * \subsection minimal_inputfile minimal.input
 * \include minimal.input
 * \section minimal_out Output
 * \include examples/minimal.out
 */
#include <faunus/faunus.h>
using namespace Faunus;                               // use Faunus namespace
typedef Geometry::Cuboid Tgeo;                        // select simulation geometry and pair potential
typedef Potential::CombinedPairPotential<Potential::Coulomb, Potential::LennardJones> Tpair;
int main() {
  atom.includefile("minimal.atoms");                  // load atom properties
  InputMap in("minimal.input");                       // open parameter file for user input
  Energy::Hamiltonian pot;                            // Hamiltonian - defines the energy field
  pot.create( Energy::Nonbonded<Tpair,Tgeo>(in) );    // add energy term for non-bonded interactions
  Space spc( pot.getGeometry() );                     // create simulation space, particles etc.
  GroupAtomic salt(spc, in);                          // group for salt particles
  Move::AtomicTranslation mv(in, pot, spc);           // particle move class
  mv.setGroup(salt);                                  // tells move class to act on salt group
  mv.move(1e5);                                       // move salt randomly 100000 times
  std::cout << spc.info() << pot.info() << mv.info(); // final information
}
