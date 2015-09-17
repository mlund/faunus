#include <faunus/faunus.h>

using namespace Faunus;

typedef Space<Geometry::Cuboidslit> Tspace;

int main() {
  cout << textio::splash();          // Spam

  InputMap mcp("slitpolymer.json");  // Open input parameter file
  EnergyDrift sys;                   // track system energy drifts
  UnitTest test(mcp);                // unit testing

  Tspace spc(mcp);                   // Simulation space (all particles and group info)
  Energy::ExternalPotential<Tspace,Potential::GouyChapman<double,true> > pot(mcp);
  pot.expot.setSurfPositionZ( &spc.geo.len_half.z() ); // Pass position of GC surface

  Move::Propagator<Tspace> mv(mcp,pot,spc);           // assembly of all MC moves
  Analysis::CombinedAnalysis analyzer(mcp,pot,spc);   // sample polymer shape etc.
  Analysis::LineDistribution<> surfmapall;            // monomer-surface histogram

  spc.load("state");                                  // Load start configuration, if any
  sys.init( Energy::systemEnergy(spc,pot,spc.p) );    // Store total system energy

  cout << spc.info() + pot.info() + textio::header("MC Simulation Begins!");

  auto pol = spc.findMolecules("polymer");

  MCLoop loop(mcp);    // handle mc loops
  while ( loop[0] ) {
    while ( loop[1] ) {
      sys += mv.move();
      analyzer.sample();

      if (slump() < 0.05) {
        for ( auto i : *pol.front() )
          surfmapall( pot.expot.surfDist( spc.p[i] ) )++;  // sample monomer distribution
      }

    } // end of micro loop

    sys.checkDrift( Energy::systemEnergy(spc,pot,spc.p) ); // re-calc system energy and detect drift
    cout << loop.timing();                                 // print timing and ETA information

  } // end of macro loop

  spc.save("state");               // save final state of simulation (positions etc)
  surfmapall.save("surfall.dat");  // save monomer-surface distribution
  FormatPQR::save("confout.pqr", spc.p);  // save PQR file

  // Perform unit tests (only for faunus integrity)
  sys.test(test);
  mv.test(test);
  analyzer.test(test);

  cout << loop.info() + sys.info() + mv.info() 
    + analyzer.info() + spc.info() + test.info();

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
