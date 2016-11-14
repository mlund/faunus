#include <faunus/faunus.h>
using namespace Faunus;
using namespace Faunus::Potential;

//typedef CombinedPairPotential<CoulombWolf,LennardJonesLB> Tpairpot; // pair potential
typedef HardSphereCap Tpairpot; // pair potential
//typedef HardSphere Tpairpot; // pair potential

typedef Geometry::Cuboid Tgeometry;   // geometry: cube w. periodic boundaries
typedef Space<Tgeometry,CapParticle> Tspace;

int main() {
  cout << textio::splash();           // show faunus banner and credits

  InputMap mcp("capparticles.json");          // open user input file
  MCLoop loop(mcp);                   // class for handling mc loops
  EnergyDrift sys;                    // class for tracking system energy drifts

  // Construct Hamiltonian and Space
  Tspace spc(mcp);

  auto pot = Energy::NonbondedVector<Tspace,Tpairpot>(mcp);

  // Markov moves and analysis
  Move::Propagator<Tspace> mv(mcp,pot,spc);
  Analysis::RadialDistribution<> rdf(0.1);      // 0.1 angstrom resolution
  Analysis::CombinedAnalysis analyzer(mcp,pot,spc);

  /*
  spc.p[0] = Point(0,0,0);
  spc.p[0].cap_center_point = Point(mcp["xyz"]["a_cap_x"],mcp["xyz"]["a_cap_y"],mcp["xyz"]["a_cap_z"]);
  spc.p[0].cap_radius = mcp["xyz"]["a_cap_radius"];
  spc.p[1] = Point(mcp["xyz"]["b_x"],mcp["xyz"]["b_y"],mcp["xyz"]["b_z"]);
  spc.p[1].cap_center_point = Point(0,0,0);
  spc.p[1].cap_radius = 0.0;
  spc.trial = spc.p;
  */
  
  for(unsigned int i = 0; i < spc.p.size(); i++)
    cout << "particle " << i << ": " << spc.p[i].is_sphere << endl;
  //return 0;
  
  spc.load("state");                               // load old config. from disk (if any)
  sys.init( Energy::systemEnergy(spc,pot,spc.p)  );// store initial total system energy

  cout << atom.info() + spc.info() + pot.info() + textio::header("MC Simulation Begins!");

  while ( loop[0] ) {  // Markov chain 
    while ( loop[1] ) {
      sys += mv.move();
      analyzer.sample();
    } // end of micro loop

    sys.checkDrift(Energy::systemEnergy(spc,pot,spc.p)); // compare energy sum with current
    cout << loop.timing();

  } // end of macro loop

  // save to disk
  FormatPQR::save("confout.pqr", spc.p); // final PQR snapshot for VMD etc.
  rdf.save("rdf.dat");                // g(r) - not normalized!
  spc.save("state");                     // final simulation state

  // perform unit tests (irrelevant for the simulation)
  UnitTest test(mcp);                    // class for unit testing
  mv.test(test);
  sys.test(test);

  // print information
  cout << loop.info() + sys.info() + mv.info() + analyzer.info() + test.info() + pot.info();
  //cout << sys.info();

  return test.numFailed();
}

/**
  @page example_capparticles Example: Cap-particles

  In this example we simulate cap-particles.

  Information about the input file can be found in `src/examples/capparticles.run`.

  capparticles.json
  =========
  @includelineno examples/capparticles.json

  capparticles.cpp
  ========
  @includelineno examples/capparticles.cpp

*/
