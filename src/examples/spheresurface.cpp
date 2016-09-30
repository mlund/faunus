#include <faunus/faunus.h>
using namespace Faunus;
using namespace Faunus::Potential;

typedef CombinedPairPotential<Coulomb,HardSphere> Tpairpot; // pair potential

typedef Geometry::SphereSurface Tgeometry;   // geometry: cube w. periodic boundaries
typedef Space<Tgeometry,PointParticle> Tspace;

int main() {
  cout << textio::splash();           // show faunus banner and credits

  InputMap mcp("spheresurface.json"); // open user input file
  MCLoop loop(mcp);                   // class for handling mc loops
  EnergyDrift sys;                    // class for tracking system energy drifts
  Tspace spc(mcp);

  auto pot = Energy::Nonbonded<Tspace,Tpairpot>(mcp);

  // Markov moves and analysis
  Move::Propagator<Tspace> mv(mcp,pot,spc);
  Analysis::RadialDistribution2D<> rdf(0.05);      // 0.05 angstrom resolution
  rdf.setRadius(mcp["system"]["spheresurface"]["radius"]); // Set radius of the 2D sphere-surface for the analysis

  spc.load("state");                               // load old config. from disk (if any)
  sys.init( Energy::systemEnergy(spc,pot,spc.p)  );// store initial total system energy

  cout << atom.info() + spc.info() + pot.info() + textio::header("MC Simulation Begins!");

  while ( loop[0] ) {  // Markov chain 
    while ( loop[1] ) {
      sys += mv.move();
      if ( slump() > -1.0 ) {
        rdf.sample( spc, atom["Pos"].id, atom["Pos"].id );
      }
    } // end of micro loop

    sys.checkDrift(Energy::systemEnergy(spc,pot,spc.p)); // compare energy sum with current
    cout << loop.timing();

  } // end of macro loop

  // save to disk
  FormatPQR::save("confout.pqr", spc.p); // final PQR snapshot for VMD etc.
  rdf.save("rdf.dat");                   // g(r) - not normalized!
  spc.save("state");                     // final simulation state

  // perform unit tests (irrelevant for the simulation)
  UnitTest test(mcp);                    // class for unit testing
  mv.test(test);
  sys.test(test);
  
  // print information
  cout << loop.info() + sys.info() + mv.info() + test.info();

  return test.numFailed();
}

/**
  @page example_bulk Example: Charges on a 2D hypersphere-surface

  The `spheresurface.cpp` program can be used to get evenly spread particles on a surface of a sphere.

  We have the following MC moves:
  - particle translation

  Information about the input file can be found in `src/examples/spheresurface.run`.

  spheresurface.json
  =========
  @includelineno examples/spheresurface.json

  spheresurface.cpp
  ========
  @includelineno examples/spheresurface.cpp

*/
