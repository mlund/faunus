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
  EnergyDrift sys;                   // class for tracking system energy drifts
  UnitTest test(mcp);                // class for unit testing (used only for Faunus integrity checks)

  // Set up energy field
  Geometry::Cuboidslit geo(mcp);     // Rectangular simulation box w. XY periodicity
  Energy::GouyChapman pot(mcp);      // Gouy-Chapman surface
  pot.setGeometry(geo);              // Pass on geometry to potential
  pot.setPosition( geo.len_half.z() ); // z position of charged surface
  Space spc( pot.getGeometry() );    // Simulation space (contains all particles and info about groups)

  // Load and add polymer to Space
  FormatAAM aam;                                      // AAM structure file I/O
  string polyfile = mcp.get<string>("polymer_file", "");
  aam.load(polyfile);                                 // load polymer structure into aam class
  Geometry::FindSpace().find(*spc.geo, spc.p, aam.particles()); // find empty spot
  GroupMolecular pol = spc.insert( aam.particles() ); // Insert into Space and return matching group
  pol.name="polymer";                                 // Give polymer arbitrary name
  spc.enroll(pol);                                    // Enroll polymer in Space

  // MC moves
  Move::TranslateRotate gmv(mcp,pot,spc);
  Move::AtomicTranslation mv(mcp,pot,spc);
  Move::CrankShaft crank(mcp,pot,spc);
  Move::Pivot pivot(mcp,pot,spc);
  Move::Reptation rep(mcp,pot,spc);

  Analysis::PolymerShape shape;                       // sample polymer shape
  Analysis::LineDistribution<> surfmapall;            // monomer-surface histogram
  spc.load("state");                                  // Load start configuration, if any
  sys.init( Energy::systemEnergy(spc,pot,spc.p) );    // Store total system energy

  cout << spc.info() << pot.info() << pol.info() << textio::header("MC Simulation Begins!");

  while ( loop.macroCnt() ) {  // Markov chain 
    while ( loop.microCnt() ) {
      int i=slp_global.rand() % 4;
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

      if (slp_global()>0.95)
        shape.sample(pol,spc);
      if (slp_global()>0.5)
        for (auto i : pol)
          surfmapall( pot.dist2surf( spc.p.at(i) ) )++;

      shape.sample(pol,spc);   // sample polymer shape - gyration radius etc.
      for (auto i : pol)       // update monomer-surface histogram
        surfmapall( pot.dist2surf( spc.p.at(i) ) )++; 

    } // end of micro loop

    sys.checkDrift( Energy::systemEnergy(spc,pot,spc.p) ); // re-calc system energy and detect drift
    cout << loop.timing();                                 // print timing and ETA information

  } // end of macro loop

  spc.save("state");               // save final state of simulation (positions etc)
  surfmapall.save("surfall.dat");  // save monomer-surface distribution
  FormatPQR().save("confout.pqr", spc.p);  // save PQR file

  // Perform unit tests (only for faunus integrity)
  gmv.test(test);
  mv.test(test);
  sys.test(test);
  pivot.test(test);
  crank.test(test);
  rep.test(test);
  shape.test(test);

  cout << loop.info() + sys.info() + gmv.info() + mv.info() + crank.info()
    + pivot.info() + rep.info() + shape.info() + spc.info() + test.info();

  return test.numFailed();
}
