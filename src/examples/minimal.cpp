#include <faunus/faunus.h>
using namespace Faunus;                        // use Faunus namespace
typedef Geometry::Cuboid Tgeo;                 // select simulation geometry and pair potential
typedef Potential::CoulombLJ Tpair;
int main() {
  ::atom.includefile("minimal.json");          // load atom properties
  InputMap in("minimal.input");                // open parameter file for user input
  Energy::Nonbonded<Tpair,Tgeo> pot(in);       // create Hamiltonian, non-bonded only
  Space spc( pot.getGeometry() );              // create simulation space, particles etc.
  Group salt;                                  // group for salt particles
  salt.addParticles(spc, in);                  // add salt according to inputfile
  Move::AtomicTranslation mv(in, pot, spc);    // particle move class
  mv.setGroup(salt);                           // tells move class to act on salt group
  mv.move(1e5);                                // move salt randomly 100000 times
  cout << spc.info() + pot.info() + mv.info(); // final information
}
/*!
 * \page example_minimal Example: Hello Monte Carlo!
 *
 * This is a minimal example of how to set up a Metropolis Monte Carlo simulation
 * with the following characteristics:
 * - Charged Lennard-Jones particles in a dielectric continuum
 * - Periodic boundaries, minimum image convention
 * - Canonical ensemble (NVT)
 * - Parameters and atom properties are read from disk
 *
 * This amounts to 15 lines of C++ code as illustated in the minimal.cpp program:
 * \includelineno minimal.cpp
 * Run the code directly from the faunus directory:
 *
 *     $ make example_minimal
 *     $ cd src/examples/
 *     $ ./minimal
 *
 * Let's walk through the code line by line:
 * - **line 1-2**
 *   - Include faunus header files and add the Faunus namespace to search path
 *
 * - **line 3**
 *   - Here we define what kind of simulation cell we wish to use - a periodic box, Faunus::Geometry::Cuboid.
 *     There are many other geometries including spheres, slits and cylinders. The geometry takes
 *     care of all distance calculations as well as boundary conditions.
 *
 * - **line 4**
 *   - This defines the pair potantial between our particles. Pair potentials can
 *     be arbitrarily mixed and Faunus::Potential::CoulombLJ is actually just a typedef
 *     for a mixture of Faunus::Potential::Coulomb and Faunus::Potential::LennardJones. To see
 *     a full list of pair potentials, check out the Faunus::Potential namespace.
 *
 * - **line 6**
 *   - Tell the atom objects (a global instance of Faunus::AtomMap) to load atomic properties
 *     from the file \ref minimal_json.
 *
 * - **line 7**
 *   - Load the user input parameter file \ref minimal_input into a Faunus::InputMap object. This simply uses keywords
 *     to get user input and is heavily used in constructors throughout Faunus.
 *
 * - **line 8**
 *   - Specify how to calculate energies in the system - i.e. the Hamiltonian. Here we only
 *     have non-bonded interactions and we need to specify the pair potential and geometry
 *     to use. Energy evaluations in faunus are done by classes derived from
 *     Faunus::Energy::Energybase.
 *     We will later see how we can construct more advanced Hamiltonians
 *     to adding Energy classes together. For a list of energy classes, see Faunus::Energy.
 *
 * - **line 9**
 *   - Faunus::Space takes care of inserting, storing and deleting particles and knows about all particle groups in the system.
 *
 * - **line 10-11**
 *   - Create a new group defining the particle range. Reads user input
 *     from the Faunus::InputMap object and automatically adds the specified particles to Space.
 *
 * - **line 12**
 *   - Instantiate a Monte Carlo move object for translating atomic particles. Moves always take care of
 *     generating a trial move, calculate the energy change, accepting/rejecting and collecting
 *     statistics. For a list of all MC moves, see Faunus::Move.
 *
 * - **line 13**
 *   - Tell the move object which particles it should move
 *
 * - **line 14**
 *   - Perform 10000 Metropolis Monte Carlo moves. Particles are randomly selected, moved, and depending on the
 *     energy change accepted or rejected. Note that the move requires access both the Hamiltonian as well as the
 *     Space. The former is used for energy evaluation, while the latter is needed to move particles. 
 *
 * - **line 15**
 *   - Print final information to standard output.
*
* If you prefer Python over C++, most of Faunus can be accessed through the `pyfaunus` module.
* The minimal.py script below is equivalent to the above C++ version:
* \includelineno minimal.py
*
* Input           {#minimalin}
* =====
* All listed files including the above C++ and python programs can be found in `src/examples/`
*
* minimal.json    {#minimal_json}
* ------------
* This input file contains properties of atoms as described in Faunus::AtomMap. For this particular
* case we set the charge ("q"), sigma ("sigma"), and displacement parameter ("dp").
* \include minimal.json
*
* minimal.input   {#minimal_input}
* -------------
* \include minimal.input
*
* Output
* ======
* This is the output generated by the `minimal.cpp` program. Note that the default is to
* use Unicode UTF-16 encoding to print mathematical symbols. If your terminal is unable to
* print this properly, Unicode output can be switched off during compilation.
* \include examples/minimal.out
*/

