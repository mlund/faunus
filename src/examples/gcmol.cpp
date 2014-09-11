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
  
  InputMap mcp("gcmol.input");     // open user input file
  MCLoop loop(mcp);                   // class for handling mc loops
  EnergyDrift sys;                    // class for tracking system energy drifts
  UnitTest test(mcp);

  // Create Space and a Hamiltonian with three terms
  Tspace spc(mcp);
  auto pot = Energy::Nonbonded<Tspace,Tpairpot>(mcp);

  string file = mcp.get<string>("polymer_file", "");


  //
  //    ADDING SPECIES, NOTE: MUST ALLOCATE GROUPS ON HEAP (use operator new)
  //
  string na("Na");
  string cl("Cl");
  string mg("Mg");

  Group *salt = new Group("salt"); // group for salt particles
  salt->setfront(spc.p.size());
  salt->addParticles(spc, na, 8);
  salt->addParticles(spc, cl, 8);
  salt->setMassCenter(spc);
  spc.enroll(*salt);

  // Add polymers
  vector<Group*> pol( mcp.get("polymer_N",0));            // vector of polymers

  for (auto* g : pol) { // load polymers
    g = new Group();
    Tspace::ParticleVector v;                            // temporary, empty particle vector
    FormatAAM::load(file,v);                             // load AAM structure into v
    Geometry::FindSpace().find(spc.geo,spc.p,v);         // find empty spot in particle vector
    *g = spc.insert(v);
    g->name="polymer";
    spc.enroll(*g);
  }

  Group *salt2 = new Group("salt2"); // group for salt particles
  salt2->setfront(spc.p.size());
  salt2->addParticles(spc, mg, 8);
  salt2->addParticles(spc, cl, 16);
  salt2->setMassCenter(spc);
  spc.enroll(*salt2);

  // Add polymers
  vector<Group*> pol2( mcp.get("polymer2_N",0));            // vector of polymers

  for (auto* g : pol2) { // load polymers
    g = new Group();
    Tspace::ParticleVector v;                            // temporary, empty particle vector
    FormatAAM::load(file,v);                             // load AAM structure into v
    Geometry::FindSpace().find(spc.geo,spc.p,v);         // find empty spot in particle vector
    *g = spc.insert(v);
    g->name="polymer2";
    spc.enroll(*g);
  }

  Group *chloride = new Group("chloride"); // group for salt particles
  chloride->setfront(spc.p.size());
  chloride->addParticles(spc, cl, 16);
  chloride->setMassCenter(spc);
  spc.enroll(*chloride);

  salt = salt2 = chloride = NULL;

  //
  //    ALL SPECIES ADDED
  //

  //
  //    ADD CONFIGURATIONS FOR POOL INSERTS
  //
  p_vec conf;
  FormatAAM::load(file,conf);
  molecule.pushConfiguration("polymer", conf);   // p_vec is one molecule large
  molecule.pushConfiguration("polymer2",conf);

  // Markov moves and analysis
  Move::GCMolecular<Tspace> gc(mcp, pot, spc);
  Move::TranslateRotate<Tspace> mv(mcp,pot,spc);

  sys.init( Energy::systemEnergy(spc,pot,spc.p)  );      // store initial total system energy

  cout << atom.info() << molecule.info() << gc.info() << spc.info() << pot.info() << textio::header("MC Simulation Begins!");

  while ( loop.macroCnt() ) {  // Markov chain
    while ( loop.microCnt() ) {
      sys+=gc.move();
    }
    sys.checkDrift(Energy::systemEnergy(spc,pot,spc.p)); // compare energy sum with current
    cout << loop.timing();
  } // end of macro loop

  // save to disk

  spc.save("state");
  FormatPQR::save("confout.pqr", spc.p);

  sys.test(test);
  gc.test(test);

  // print information
  cout << loop.info() << sys.info() << gc.info() << spc.info() << test.info();

  //
  // clean allocated memory
  //
  spc.freeGroups();

  return test.numFailed();
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


