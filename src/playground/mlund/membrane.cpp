#include <faunus/faunus.h>
#include <faunus/membrane.h>

using namespace Faunus;
using namespace Faunus::Potential;

typedef Geometry::Cuboid Tgeometry;   // specify geometry - here cube w. periodic boundaries
//typedef CosbinedPairPotential<WeeksChandlerAndersen, DebyeHuckel> Tdhwca;
typedef CosAttractCombi<WeeksChandlerAndersen> Tpairpot;

int main() {
  cout << textio::splash();           // show faunus banner and credits

  InputMap mcp("membrane.input");     // open user input file
  MCLoop loop(mcp);                   // class for handling mc loops
  FormatPQR pqr;                      // PQR structure file I/O
  FormatAAM aam;
  FormatXTC xtc(1000);                // XTC gromacs trajectory format
  EnergyDrift sys;                    // class for tracking system energy drifts
  UnitTest test(mcp);                 // class for unit testing

  // Energy functions and space
  Energy::Hamiltonian pot;
  auto nonbonded = pot.create( Energy::Nonbonded<Tpairpot,Tgeometry>(mcp) );
  auto bonded    = pot.create( Energy::Bonded() );
  Space spc( pot.getGeometry() );

  // Markov moves and analysis
  Move::TranslateRotate gmv(mcp,pot,spc);
  Move::AtomicTranslation mv(mcp, pot, spc);
  Move::Pivot pivot(mcp, pot, spc);
  Move::Reptation rep(mcp, pot, spc);
  Move::CrankShaft crank(mcp, pot, spc);
  Analysis::PolymerShape shape;
  Analysis::RadialDistribution<> rdf(0.2);

  // Load membrane
  DesernoMembrane<Tgeometry> mem(mcp,pot,spc, nonbonded->pairpot.first, nonbonded->pairpot.second);

  // Load peptide
  string polyfile = mcp.get<string>("polymer_file", "");
  double req    = mcp.get<double>("polymer_eqdist", 0);
  double k      = mcp.get<double>("polymer_forceconst", 0);
  aam.load(polyfile);                                 // load polymer structure into aam class
  Geometry::FindSpace f;                              // class for finding empty space in container
  f.find(*spc.geo, spc.p, aam.particles());           // find empty spot
  GroupMolecular pol = spc.insert( aam.particles() ); // Insert particles into Space and return matching group
  pol.name="peptide";                                 // Give the polymer an arbitrary name
  spc.enroll(pol);                                    // All groups need to be enrolled in the Space
  //for (int i=pol.front(); i<pol.back(); i++)
  //  bonded->add(i, i+1, Potential::Harmonic(k,req));   // add bonds

  spc.load("state");                                     // load old config. from disk (if any)
  pqr.save("initial.pqr", spc.p);

  sys.init( Energy::systemEnergy(spc,pot,spc.p)  );      // store initial total system energy

  cout << atom.info() << spc.info() << pot.info() << textio::header("MC Simulation Begins!");

  while ( loop.macroCnt() ) {  // Markov chain 
    while ( loop.microCnt() ) {
      int k=mem.lipids.sizeMol(); //number of lipids
      int i=slp_global.rand() % 3;
      Group g;
      switch (i) {
        case 0:
          mv.setGroup(mem.lipids);
          sys+=mv.move( mem.lipids.size() ); // translate lipid monomers
          break;
        case 1:
          while (k-->0) {
            i=mem.lipids.randomMol(); // pick random lipid molecule
            g=mem.lipids[i];
            g.name="subgrouplipid";
            g.setMassCenter(spc);       // mass center needed for rotation
            gmv.setGroup(g);            // tell what to move
            sys+=gmv.move();            // translate/rotate polymers
          }
          break;
        case 2:
          gmv.setGroup(pol);
          gmv.move();
      }
      if ( slp_global.runtest(0.1) )
        xtc.save("traj.xtc", spc);

    } // end of micro loop

    sys.checkDrift(Energy::systemEnergy(spc,pot,spc.p)); // compare energy sum with current

    // save to disk
    rdf.save("rdf.dat");
    pqr.save("confout.pqr", spc.p);
    spc.save("state");

    cout << loop.timing();

  } // end of macro loop

  // perform unit tests
  gmv.test(test);
  mv.test(test);
  sys.test(test);

  // print information
  cout << loop.info() << sys.info() << mv.info() << gmv.info() << rep.info()
    << pivot.info() << shape.info() << test.info() << spc.info();

  return test.numFailed();
}
