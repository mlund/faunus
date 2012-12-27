#include <faunus/faunus.h>
using namespace Faunus;                               // use Faunus namespace
typedef Geometry::Cuboid Tgeo;                        // select simulation geometry and pair potential
typedef Potential::CoulombLJ Tpair;
int main() {
  atom.includefile("minimal.json");                   // load atom properties
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
/*!
 * \page example_minimal Example: Hello Monte Carlo!
 * This is a minimal example of how to set up a Metropolis Monte Carlo of atomic
 * particles with the following characteristics:
 * - Charged Lennard-Jones particles in a dielectric continuum
 * - Periodic boundaries, minimum image convention
 * - Canonical ensemble (NVT)
 *
 * This amounts to 16 lines of C++ code and 12 lines of python
 * code as shown below.
 * Let us now walk through the C++ code line by line:
 * - <b>line 1-2</b>
 *   - Include faunus header files and add the Faunus namespace to search path
 * - <b>line 3</b>
 *   - Here we define what kind of simulation cell we wish to use - a periodic box, Faunus::Geometry::Cuboid.
 *     There are many other geometries including spheres, slits and cylinders. The geometry takes
 *     care of all distance calculations as well as boundary conditions.
 * - <b>line 4</b>
 *   - This defines the pair potantial between our particles. Pair potentials can
 *     be arbitrarily mixed and Faunus::Potential::CoulombLJ is actually just a typedef
 *     for a mixture of Faunus::Potential::Coulomb and Faunus::Potential::LennardJones. To see
 *     a full list of pair potentials, check out the Faunus::Potential namespace.
 * - <b>line 6</b>
 *   - Tell the atom objects (a global instance of Faunus::AtomTypes) to load atomic properties
 *     from the file \c "minimal.json".
 * - <b>line 7</b>
 *   - Load the user input parameter file \c "minimal.input" into a Faunus::InputMap object. This simply uses keywords
 *     to get user input and is heavily used in constructors throughout Faunus.
 * - <b>line 8</b>
 *   - Create a Hamiltonian that can sum interaction energies from an arbitrary number of energy classes.
 *     Energy classes evaluate specific contributions to the total energy -- for examples non-bonded interactions,
 *     bonded interactions, external potentials etc. For a list of all energy classes, see Faunus::Energy.
 * - <b>line 9</b>
 *   - Add nonbonded energy class to the Hamiltonian. The Faunus::Energy::Nonbonded template is contructed with the
 *     geometry and pair potential as these are needed to evaluate energies. Although we could add many
 *     energy classes to the Hamiltonian, we here suffice with merely one.
 * - <b>line 10</b>
 *   - Faunus::Space takes care of inserting, storing and deleting particles and knows about all particle groups in the system.
 * - <b>line 11</b>
 *   - Create a new group defining the particle range. The constructor of Faunus::GroupAtomic reads user input
 *     from the Faunus::InputMap object and automatically adds the specified particles to Space.
 * - <b>line 12</b>
 *   - Instantiate a Monte Carlo move object for translating atomic particles. Moves always take care of
 *     generating a trial move, calculate the energy change, accepting/rejecting and collecting
 *     statistics. For a list of all MC moves, see Faunus::Move.
 * - <b>line 13</b>
 *   - Tell the move object which particles it should move
 * - <b>line 14</b>
 *   - Perform 10000 Metropolis Monte Carlo moves. Particles are randomly selected, moved, and depending on the
 *     energy change accepted or rejected. Note that the move requires access both the Hamiltonian as well as the
 *     Space. The former is used for energy evaluation, while the latter is needed to move particles. 
 * - <b>line 15</b>
 *   - Print final information to standard output.
 *
 * \section minimal_source Source files
 * All of the files listed below can be found in the \c src/examples/ directory.
 * \subsection minimal_cpp minimal.cpp
 * \includelineno minimal.cpp
 * \subsection minimal_py minimal.py
 * \includelineno minimal.py
 * \section minimal_input Input
 * \subsection minimal_json minimal.json
 * This input file contains properties of atoms as described in Faunus::AtomTypes. For this particular
 * case we set the charge ("q"), sigma ("sigma"), and displacement parameter ("dp").
 * \include minimal.json
 * \subsection minimal_inputfile minimal.input
 * \include minimal.input
 * \section minimal_out Output
 * \include examples/minimal.out
 */

