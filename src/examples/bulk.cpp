#include <faunus/faunus.h>
#include <faunus/ewald.h>
using namespace Faunus;
using namespace Faunus::Potential;

#if defined(COULOMB)
typedef CombinedPairPotential<Coulomb,LennardJonesTrunkShift> Tpairpot; // pair potential
#elif defined(DEBYEHUCKEL)
typedef CombinedPairPotential<DebyeHuckelDenton,LennardJonesTrunkShift> Tpairpot; // pair potential
#else
typedef CombinedPairPotential<CoulombWolf,LennardJonesLB> Tpairpot; // pair potential
#endif

typedef Geometry::Cuboid Tgeometry;   // geometry: cube w. periodic boundaries
typedef Space<Tgeometry,PointParticle> Tspace;

int main() {
  cout << textio::splash();           // show faunus banner and credits

  InputMap mcp("bulk.input");         // open user input file
  MCLoop loop(mcp);                   // class for handling mc loops
  EnergyDrift sys;                    // class for tracking system energy drifts

  // Construct Hamiltonian and Space
  Tspace spc(mcp);
  auto pot = Energy::Nonbonded<Tspace,Tpairpot>(mcp)
    + Energy::ExternalPressure<Tspace>(mcp);

  pot.setSpace(spc);

  // Markov moves and analysis
  Move::Isobaric<Tspace> iso(mcp,pot,spc);
  Move::AtomicTranslation<Tspace> mv(mcp,pot,spc);
  Analysis::RadialDistribution<> rdf_ab(0.1);      // 0.1 angstrom resolution
  Analysis::VirialPressure virial;
  Average<double> pm;

  // Add salt
  Group salt;
  salt.addParticles(spc, mcp);
  mv.setGroup(salt);

  spc.load("state");                               // load old config. from disk (if any)
  sys.init( Energy::systemEnergy(spc,pot,spc.p)  );// store initial total system energy

  cout << atom.info() + spc.info() + pot.info()
    + textio::header("MC Simulation Begins!");

#ifdef DEBYEHUCKEL
  cout << pot.first.pairpot.first.info(spc.p,spc.geo.getVolume()) << endl;
#endif

  while ( loop.macroCnt() ) {  // Markov chain 
    while ( loop.microCnt() ) {
      if (slp_global() < 0.5)
        sys+=mv.move( salt.size() );  // translate salt
      else 
        sys+=iso.move();              // isobaric volume move

      if (slp_global() < 0.05) {
        rdf_ab.sample(spc,salt,atom["Na"].id,atom["Cl"].id);
        virial.sample(spc,pot);
        //pm+=pot.first.pairpot.first.virial(spc.p,spc.geo);
      }
    } // end of micro loop

    sys.checkDrift(Energy::systemEnergy(spc,pot,spc.p)); // compare energy sum with current
    cout << loop.timing();

  } // end of macro loop

  // save to disk
  FormatPQR::save("confout.pqr", spc.p); // final PQR snapshot for VMD etc.
  rdf_ab.save("rdf.dat");                // g(r) - not normalized!
  spc.save("state");                     // final simulation state

  // perform unit tests (irrelevant for the simulation)
  UnitTest test(mcp);                    // class for unit testing
  iso.test(test);
  mv.test(test);
  sys.test(test);

  // print information
  cout << loop.info() + sys.info() + mv.info() + iso.info()
    + virial.info() + test.info();

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
  The `bulk.cpp` program can be used to simulate any atomic mixtures and the
  dielectric constant may also be varied (it is unity in this example).

  We have the following MC moves:
  - salt translation
  - isotropic volume move (NPT ensemble)

  Information about the input file can be found in `src/examples/bulk.run`.

  ![Na-Cl distribution function with various electrostatic potentials.](wolf.png)

  bulk.cpp
  ========
  @includelineno examples/bulk.cpp

*/
