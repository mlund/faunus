#include <faunus/faunus.h>

using namespace Faunus;

int main() {
  cout << textio::splash();          // Spam

  InputMap mcp("slitpolymer.input"); // Open input parameter file
  MCLoop loop(mcp);                  // handle mc loops
  EnergyDrift sys;                   // track system energy drifts
  UnitTest test(mcp);                // unit testing
  Geometry::Cuboidslit geo(mcp);     // rectangular simulation box w. XY periodicity

  Energy::ExternalPotential< Potential::GouyChapman<> > pot(mcp);
  pot.setGeometry(geo);              // Pass on geometry to potential
  pot.expot.setSurfPositionZ( &geo.len_half.z() ); // Pass position of GC surface
  Space spc( pot.getGeometry() );    // Simulation space (all particles and group info)

  // Load and add polymer to Space
  FormatAAM aam;                                      // AAM structure file I/O
  string polyfile = mcp.get<string>("polymer_file", "");
  aam.load(polyfile);                                 // Load polymer structure into aam class
  Geometry::FindSpace().find(*spc.geo, spc.p, aam.particles()); // find empty spot
  GroupMolecular pol = spc.insert( aam.particles() ); // Insert into Space and return matching group
  pol.name="polymer";                                 // Give polymer arbitrary name
  spc.enroll(pol);                                    // Enroll polymer in Space

  // MC moves
  Move::TranslateRotate gmv(mcp,pot,spc);
  Move::CrankShaft crank(mcp,pot,spc);
  Move::Pivot pivot(mcp,pot,spc);
  Move::Reptation rep(mcp,pot,spc);

  Analysis::PolymerShape shape;                       // sample polymer shape
  Analysis::LineDistribution<> surfmapall;            // monomer-surface histogram
  spc.load("state");                                  // Load start configuration, if any
  sys.init( Energy::systemEnergy(spc,pot,spc.p) );    // Store total system energy

  cout << spc.info() + pot.info() + pol.info() + textio::header("MC Simulation Begins!");

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

      double rnd = slp_global(); // [0:1[
      if (rnd<0.05)
        shape.sample(pol,spc);   // sample polymer shape - gyration radius etc.
      if (rnd<0.05)
        for (auto i : pol)
          surfmapall( pot.expot.surfDist( spc.p[i] ) )++;  // sample monomer distribution

    } // end of micro loop

    sys.checkDrift( Energy::systemEnergy(spc,pot,spc.p) ); // re-calc system energy and detect drift
    cout << loop.timing();                                 // print timing and ETA information

  } // end of macro loop

  spc.save("state");               // save final state of simulation (positions etc)
  surfmapall.save("surfall.dat");  // save monomer-surface distribution
  FormatPQR().save("confout.pqr", spc.p);  // save PQR file

  // Perform unit tests (only for faunus integrity)
  gmv.test(test);
  sys.test(test);
  pivot.test(test);
  crank.test(test);
  rep.test(test);
  shape.test(test);

  cout << loop.info() + sys.info() + gmv.info() + crank.info()
    + pivot.info() + rep.info() + shape.info() + spc.info() + test.info();

  return test.numFailed();
}

/**
  @page example_slitpolymer Example: Polymer/Gouy-Chapman surface

  We simulate a single polyelectrolyte in a rectangular slit where one wall
  is charged. The following interactions are included,

  - monomer-surface: Gouy-Chapman electrostatics
  - monomer-monomer: stiff bonds - i.e. ideal chain

  This is naturally an unphysical system but can be compared with *exact*
  polymer DFT theory -- see `slitpolymer.agr` xmgrace plot.

  slitpolymer.cpp
  ===============

  @includelineno examples/slitpolymer.cpp

*/
