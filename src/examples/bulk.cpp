#include <faunus/faunus.h>
#include <faunus/ewald.h>
//#define EWALD
using namespace Faunus;
using namespace Faunus::Potential;

#if defined(COULOMB)
typedef CombinedPairPotential<CoulombGalore,LennardJonesTrunkShift> Tpairpot; // pair potential
#elif defined(DEBYEHUCKEL)
typedef CombinedPairPotential<DebyeHuckelDenton,LennardJonesTrunkShift> Tpairpot; // pair potential
#elif defined(EWALD)
typedef LennardJonesLB Tpairpot;
#else
typedef CombinedPairPotential<CoulombGalore,LennardJonesLB> Tpairpot; // pair potential
typedef CutShift<Tpairpot,false> TpairpotCut;
#endif

typedef Geometry::Cuboid Tgeometry;   // geometry: cube w. periodic boundaries
typedef Space<Tgeometry,PointParticle> Tspace;

int main() {
  cout << textio::splash();           // show faunus banner and credits

  Tmjson mcp = openjson("bulk.json"); // open JSON input file
  MCLoop loop(mcp);                   // class for handling mc loops

  Tspace spc(mcp);                    // simulation space

#ifdef EWALD
  auto pot = Energy::NonbondedEwald<Tspace,Tpairpot>(mcp)
    + Energy::ExternalPressure<Tspace>(mcp);
#else
  auto pot = Energy::NonbondedCutg2g<Tspace,TpairpotCut>(mcp)
    + Energy::ExternalPressure<Tspace>(mcp);
#endif

  spc.load("state");                  // load old config. from disk (if any)

  Analysis::CombinedAnalysis analyzer(mcp,pot,spc);
  Move::Propagator<Tspace> mv(mcp,pot,spc);

  cout << atom.info() + spc.info() + pot.info()
    + textio::header("MC Simulation Begins!");

#ifdef DEBYEHUCKEL
  cout << pot.first.pairpot.first.info(spc.p,spc.geo.getVolume()) << endl;
#endif

  while ( loop[0] ) {  // Markov chain 
    while ( loop[1] ) {
      mv.move();
      analyzer.sample();
    } // end of micro loop

    cout << loop.timing();

  } // end of macro loop

  // perform unit tests (irrelevant for the simulation)
  UnitTest test(mcp);                    // class for unit testing
  mv.test(test);
  analyzer.test(test);

  // print information
  cout << loop.info() + mv.info() + analyzer.info() + test.info();

  return test.numFailed();
}

/**
  @page example_bulk Example: Melted NaCl

  In this example we simulate melted NaCl in the NVT or NPT ensemble. We use a
  Lennard-Jones potential combined with a shifted Coulombic potential according to
  [Wolf/Yonezawa](<http://dx.doi.org/10/j97>).
  This gives essentially identical results to the more elaborate Particle
  Mesh Ewald method - see figure below. In contrast, using the simple minimum image
  approach with a cubic cutoff, the system freezes.
  The `bulk.cpp` program can be used to simulate any atomic (or molecular) mixtures
  and the dielectric constant may also be varied (it is unity in this example).

  We have the following MC moves:
  - salt translation
  - isotropic volume move (NPT ensemble)

  Information about the input file can be found in `src/examples/bulk.py`.

  ![Na-Cl distribution function with various electrostatic potentials.](wolf.png)

  bulk.cpp
  ========
  @includelineno examples/bulk.cpp
*/
