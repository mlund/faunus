#include <faunus/faunus.h>

using namespace Faunus;

typedef Space<Geometry::Cuboidslit> Tspace;

int main() {
  cout << textio::splash();          // Spam

  InputMap mcp("slitpolymer.input"); // Open input parameter file
  MCLoop loop(mcp);                  // handle mc loops
  EnergyDrift sys;                   // track system energy drifts
  UnitTest test(mcp);                // unit testing

  Tspace spc(mcp);                   // Simulation space (all particles and group info)
  Energy::ExternalPotential<Tspace,Potential::GouyChapman<double,true> > pot(mcp);
  pot.expot.setSurfPositionZ( &spc.geo.len_half.z() ); // Pass position of GC surface

  // Load and add polymer to Space
  string file = mcp.get<string>("polymer_file", "");
  Tspace::ParticleVector v;                   // temporary, empty particle vector
  FormatAAM::load(file,v);                    // load AAM structure into v
  Geometry::FindSpace().find(spc.geo,spc.p,v);// find empty spot in particle vector
  Group pol = spc.insert(v);                  // Insert into Space and return matching group
  pol.name="polymer";                         // Give polymer arbitrary name
  spc.enroll(pol);                            // Enroll polymer in Space

  // MC moves
  Move::TranslateRotate<Tspace> gmv(mcp,pot,spc);
  Move::CrankShaft<Tspace> crank(mcp,pot,spc);
  Move::Pivot<Tspace> pivot(mcp,pot,spc);
  Move::Reptation<Tspace> rep(mcp,pot,spc);

  Analysis::PolymerShape shape;                       // sample polymer shape
  Analysis::LineDistribution<> surfmapall;            // monomer-surface histogram
  spc.load("state");                                  // Load start configuration, if any
  sys.init( Energy::systemEnergy(spc,pot,spc.p) );    // Store total system energy

  cout << spc.info() + pot.info() + pol.info() + textio::header("MC Simulation Begins!");

  while ( loop[0] ) {  // Markov chain 
    while ( loop[1] ) {
      int i= slump.range(0,3);
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

      double rnd = slump(); // [0:1[
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
  FormatPQR::save("confout.pqr", spc.p);  // save PQR file

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
