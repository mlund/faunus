#include <faunus/faunus.h>

using namespace Faunus;                   // use Faunus namespace
typedef Space<Geometry::Cuboid> Tspace;   // Type of simulation space
typedef Potential::LennardJonesLB Tpair;  // and pair potential

int main()
{
    InputMap in("minimal.json");            // open parameter file for user input
    Tspace spc(in);                         // Simulation space, particles etc.
    Energy::Nonbonded<Tspace, Tpair> pot(in);// Hamiltonian, non-bonded only
    Move::Propagator<Tspace> mv(in, pot, spc);// Monte Carlo move class
    for ( int i = 0; i < 1e4; i++ )
        mv.move();                            // move salt randomly N times
    cout << spc.info() + pot.info() + mv.info(); // final information
}
/**
  @page example_minimal Example: Hello Monte Carlo!
 
  This is a minimal example of how to set up a Metropolis Monte Carlo simulation
  with the following characteristics:
  - Lennard-Jones particles
  - Periodic boundaries, minimum image convention
  - Canonical ensemble (NVT)
  - Parameters and atom properties are read from disk
 
  This amounts to 14 lines of C++ code as illustrated in the minimal.cpp
  program:
 
  @includelineno minimal.cpp
 
  Run the code directly from the faunus directory:
 
      $ make example_minimal
      $ cd src/examples/
      $ ./minimal
 
  Let's walk through the code line by line:
  - **line 1-2**
    - Include faunus header files and add the Faunus namespace to search path
 
  - **line 3**
    - Here we define what kind of simulation cell we wish to use
      a periodic box, `Faunus::Geometry::Cuboid`.
      There are many other geometries including spheres, slits and
      cylinders. The geometry takes
      care of all distance calculations as well as boundary conditions.
 
  - **line 4**
    - This defines the pair potantial between our particles.
      Pair potentials can
      be arbitrarily mixed and `Faunus::Potential::CoulombLJ` is actually
      just a typedef
      for a mixture of `Faunus::Potential::Coulomb` and
      `Faunus::Potential::LennardJones`. To see
      a full list of pair potentials, check out the `Faunus::Potential`
      namespace.
 
  - **line 5**
    - This defines the simulation space, which includes all properties as
      well as geometry information
 
  - **line 7**
    - Load the user input parameters into a `Faunus::InputMap` object.
    This simply uses the JSON format
      to get user input and is heavily used in constructors throughout Faunus.
 
  - **line 8**
    - Specify how to calculate energies in the system - i.e. the Hamiltonian.
      Here we only have non-bonded interactions and we need to specify
      the pair potential and space type
      to use. Energy evaluations in faunus are done by classes derived from
      `Faunus::Energy::Energybase`.
      We will later see how we can construct more advanced Hamiltonians
      to adding Energy classes together. For a list of energy classes, see Faunus::Energy.
 
  - **line 9**
    - Space takes care of inserting, storing and deleting particles and knows about
      all particle groups in the system.
 
  - **line 10**
    - Instantiate a Monte Carlo move object for translating atomic particles.
      Moves always take care of generating a trial move, calculate the energy change,
      accepting/rejecting and collecting statistics. For a list of all MC moves,
      see Faunus::Move.
      The `Move::Propagator` class is a general move class that is controlled
      from the json input file. In this example we add simple translational
      move.
 
  - **line 11-12**
    - Perform N-particle Metropolis Monte Carlo moves.
      Particles are randomly selected, moved, and depending on the
      energy change accepted or rejected. Note that the move requires access to
      both the Hamiltonian as well as the
      Space. The former is used for energy evaluation, while the latter is needed
      to move particles. 
 
  - **line 13**
    - Print final information to standard output.

  Exercise
  --------
  1. Discuss why the acceptance (acc) and mean square displacement (msq)
     differ between particle A and B?

  2. Change the translational displacement parameter (dp) and observe what
     happens to the acceptance (acc) and the particle mean square displacement (msq).
     Plot (msq) vs. (dp) and explain why it is desirable to maximize
     this function.

  3. Modify the code so that a final PQR snapshot is saved at the of
     the simulation. Usuful for visualization in i.e. VMD.
     Hint: Look for `FormatPQR` in the documentation.

  4. Modify the input so that the system freezes. There are several ways to do this;
     discuss these.
 
  Input           {#minimalin}
  =====
  All listed files including the above C++ and python programs can be found in
  `src/examples/`
 
  minimal.json    {#minimal_json}
  ------------
  This input file contains properties of atoms as described in `Faunus::AtomData`;
  ditto for molecules, see `Faunus::MoleculeData`. Monte Carlo moves are specified
  in the `moves` -- see i.e. `Faunus::Move::Propagator` for a list of possible moves.

  @include minimal.json
 
  Output
  ======
  This is the output generated by the `minimal.cpp` program. Note that the default is to
  use Unicode UTF-16 encoding to print mathematical symbols. If your terminal is unable to
  print this properly, Unicode output can be switched off during compilation.
 
  @include examples/minimal.out
*/

