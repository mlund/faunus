/*
 * We simulate a single polyelectrolyte in a rectangular slit where one wall
 * is charged. The following interactions are included,
 * 1) monomer-surface: Gouy-Chapman electrostatics
 * 2) monomer-monomer: stiff bonds - i.e. ideal chain
 *
 * This is of course an unphysical system but can be calculated exactly
 * using polymer DFT theory -- see slitpolymer.agr xmgrace plot.
 */

#include <faunus/faunus.h>

using namespace Faunus;

int main() {
  cout << textio::splash();          // show spam...

  InputMap mcp("slitpolymer.input"); // open input parameter file
  MCLoop loop(mcp);                  // class for handling mc loops
  FormatPQR pqr;                     // PQR structure file I/O
  FormatAAM aam;                     // AAM structure file I/O
  EnergyDrift sys;                   // class for tracking system energy drifts
  UnitTest test(mcp);                // class for unit testing (used only for Faunus integrity checks)

  // Set up energy field
  Geometry::Cuboidslit geo(mcp);     // Rectangular slit simulation container w. XY periodicity
  Energy::GouyChapman pot(mcp);      // Gouy-Chapman electrostatics from charged surface
  pot.setGeometry(geo);              // Pass on geometry to potential
  pot.setPosition( geo.len_half.z ); // z position of charged surface
  Space spc( pot.getGeometry() );    // Simulation space (contains all particles and info about groups)

  // Load and add polymer to Space
  string polyfile = mcp.get<string>("polymer_file", "");
  aam.load(polyfile);                                 // load polymer structure into aam class
  Geometry::FindSpace f;                              // class for finding empty space in container
  f.find(*spc.geo, spc.p, aam.particles());           // find empty spot
  GroupMolecular pol = spc.insert( aam.particles() ); // Insert particles into Space and return matching group
  pol.name="polymer";                                 // Give the polymer an arbitrary name
  spc.enroll(pol);                                    // All groups need to be enrolled in the Space

  // MC moves
  Move::TranslateRotate gmv(mcp,pot,spc);
  Move::AtomicTranslation mv(mcp,pot,spc);
  Move::CrankShaft crank(mcp,pot,spc);
  Move::Pivot pivot(mcp,pot,spc);
  Move::Reptation rep(mcp,pot,spc);

  Analysis::PolymerShape shape;                       // class for sampling the polymer shape
  Analysis::LineDistribution<float,int> surfmapall;   // histogram for monomer-surface distribution

  spc.load("state");                                  // Load start configuration, if any
  sys.init( Energy::systemEnergy(spc,pot,spc.p) );    // Calculate initial, total system energy

  cout << spc.info() << pot.info() << pol.info() << textio::header("MC Simulation Begins!");

  while ( loop.macroCnt() ) {  // Markov chain 
    while ( loop.microCnt() ) {
      int i=rand() % 4;
      switch (i) {
        case 0: // translate and rotate polymer
          gmv.setGroup(pol);
          sys+=gmv.move(); 
          break;
        case 1: // pivot
          pivot.setGroup(pol);
          sys+=pivot.move();
          break;
        case 2: // crankshaft
          crank.setGroup(pol);
          sys+=crank.move();
          break;
        case 3: // reptation
          rep.setGroup(pol);
          sys+=rep.move();
          break;
      }

      shape.sample(pol,spc);   // sample polymer shape - gyration radius etc.
      for (auto i : pol)       // update monomer-surface histogram
        surfmapall( pot.dist2surf( spc.p.at(i) ) )++; 

    } // end of micro loop

    sys.checkDrift( Energy::systemEnergy(spc,pot,spc.p) ); // re-calc system energy and detect drift
    cout << loop.timing();                                 // print timing and ETA information

  } // end of macro loop

  pqr.save("confout.pqr", spc.p);  // save PQR file
  spc.save("state");               // save final state of simulation (positions etc)
  surfmapall.save("surfall.dat");  // save monomer-surface distribution

  // Perform unit tests (only for faunus integrity)
  gmv.test(test);
  mv.test(test);
  sys.test(test);
  pivot.test(test);
  crank.test(test);
  rep.test(test);
  shape.test(test);

  cout << loop.info() << sys.info() << gmv.info() << mv.info() << crank.info()
    << pivot.info() << rep.info() << shape.info() << spc.info() << test.info();

  return test.numFailed();
}
