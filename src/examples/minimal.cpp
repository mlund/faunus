#include <faunus/faunus.h>
using namespace Faunus;                   // use Faunus namespace
typedef Space<Geometry::Cuboid> Tspace;   // Type of simulation space
typedef Potential::CoulombLJ Tpair;       // and pair potential

int main() {
  ::atom.includefile("minimal.json");     // load atom properties
  InputMap in("minimal.input");           // open parameter file for user input
  Energy::Nonbonded<Tspace,Tpair> pot(in);// Hamiltonian, non-bonded only
  Tspace spc(in);                         // Simulation space, particles etc.
  Group salt;                             // Group for salt particles
  salt.addParticles(spc,in);              // Add according to user input
  Move::AtomicTranslation<Tspace> mv(in,pot,spc);// particle move class
  mv.setGroup(salt);                      // move class acts on salt group
  mv.move(1e5);                           // move salt randomly 100000 times
  cout << spc.info() + pot.info() + mv.info(); // final information
}
/**
 * @page example_minimal Example: Hello Monte Carlo!
 *
 * This is a minimal example of how to set up a Metropolis Monte Carlo simulation
 * with the following characteristics:
 * - Charged Lennard-Jones particles in a dielectric continuum
 * - Periodic boundaries, minimum image convention
 * - Canonical ensemble (NVT)
 * - Parameters and atom properties are read from disk
 *
 * This amounts to 16 lines of C++ code as illustated in the minimal.cpp
 * program:
 * @includelineno minimal.cpp
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
 *   - Here we define what kind of simulation cell we wish to use
 *     a periodic box, Faunus::Geometry::Cuboid.
 *     There are many other geometries including spheres, slits and
 *     cylinders. The geometry takes
 *     care of all distance calculations as well as boundary conditions.
 *
 * - **line 4**
 *   - This defines the pair potantial between our particles.
 *     Pair potentials can
 *     be arbitrarily mixed and `Faunus::Potential::CoulombLJ` is actually
 *     just a typedef
 *     for a mixture of `Faunus::Potential::Coulomb` and
 *     `Faunus::Potential::LennardJones`. To see
 *     a full list of pair potentials, check out the `Faunus::Potential`
 *     namespace.
 *
 * - **line 5**
 *   - This defines the simulation space, which includes all properties as
 *     well as geometry information
 *
 * - **line 7**
 *   - Tell the atom objects (a global instance of `Faunus::AtomMap`)
 *     to load atomic properties from the file \ref minimal_json.
 *
 * - **line 8**
 *   - Load the user input parameter file \ref minimal_input into a
 *     `Faunus::InputMap object`. This simply uses keywords
 *     to get user input and is heavily used in constructors throughout Faunus.
 *
 * - **line 9**
 *   - Specify how to calculate energies in the system - i.e. the Hamiltonian.
 *     Here we only have non-bonded interactions and we need to specify
 *     the pair potential and space type
 *     to use. Energy evaluations in faunus are done by classes derived from
 *     `Faunus::Energy::Energybase`.
 *     We will later see how we can construct more advanced Hamiltonians
 *     to adding Energy classes together. For a list of energy classes, see Faunus::Energy.
 *
 * - **line 10**
 *   - Space takes care of inserting, storing and deleting particles and knows about
 *     all particle groups in the system.
 *
 * - **line 11-12**
 *   - Create a new group defining the particle range. Reads user input
 *     from the Faunus::InputMap object and automatically adds the specified particles to Space.
 *
 * - **line 13**
 *   - Instantiate a Monte Carlo move object for translating atomic particles. Moves always take care of
 *     generating a trial move, calculate the energy change, accepting/rejecting and collecting
 *     statistics. For a list of all MC moves, see Faunus::Move.
 *
 * - **line 14**
 *   - Tell the move object which particles it should move
 *
 * - **line 15**
 *   - Perform 10000 Metropolis Monte Carlo moves. Particles are randomly selected, moved, and depending on the
 *     energy change accepted or rejected. Note that the move requires access both the Hamiltonian as well as the
 *     Space. The former is used for energy evaluation, while the latter is needed to move particles. 
 *
 * - **line 16**
 *   - Print final information to standard output.
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

