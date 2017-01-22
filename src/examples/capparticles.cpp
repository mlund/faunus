#include <faunus/faunus.h>
using namespace Faunus;
using namespace Faunus::Potential;

//typedef CombinedPairPotential<CoulombWolf,LennardJonesLB> Tpairpot; // pair potential
typedef HardSphereCap Tpairpot; // pair potential
//typedef CombinedPairPotential<HardSphereCap,CoulombCap> Tpairpot; // pair potential

typedef Geometry::Cuboid Tgeometry;   // geometry: cube w. periodic boundaries
typedef Space<Tgeometry,CapParticle> Tspace;

int main() {
  //cout << textio::splash();           // show faunus banner and credits

  InputMap mcp("capparticles.json");          // open user input file
  MCLoop loop(mcp);                   // class for handling mc loops
  EnergyDrift sys;                    // class for tracking system energy drifts

  // Construct Hamiltonian and Space
  Tspace spc(mcp);

  auto pot = Energy::NonbondedVector<Tspace,Tpairpot>(mcp);// + Energy::ExternalPressure<Tspace>(mcp);;

  // Markov moves and analysis
  Move::Propagator<Tspace> mv(mcp,pot,spc);
  Analysis::RadialDistribution<> rdf(0.1);      // 0.1 angstrom resolution
  Analysis::CombinedAnalysis analyzer(mcp,pot,spc);
  Analysis::CapAnalysis capa;
  FormatXTC xtc(1000);
  
  double scaled_center_point = mcp["xyz"]["scaled_center_point"];
  double scaled_radius = mcp["xyz"]["scaled_radius"];
  spc.p[0] = Point(0,0,0);
  spc.p[0].cap_center_point = Point(spc.p[0].radius*scaled_center_point,0,0);
  spc.p[0].cap_center = spc.p[0].radius*scaled_center_point;
  spc.p[0].cap_radius = spc.p[0].radius*scaled_radius;
  spc.p[0].is_sphere = false;
  spc.p[0].update();
  spc.p[1] = Point(0,0,mcp["xyz"]["z"]);
  spc.p[1].cap_center_point = Point(0,0,0);
  spc.p[1].cap_center = 0.0;
  spc.p[1].cap_radius = 0.0;
  spc.p[1].is_sphere = true;
  spc.p[1].update();
  spc.trial = spc.p;
  double minY = mcp["xyz"]["minY"];
  double maxY = mcp["xyz"]["maxY"];
  double minX = mcp["xyz"]["minX"];
  double maxX = mcp["xyz"]["maxX"];
  double step = mcp["xyz"]["step"];
  int cnt;
  for(double y = minY; y <= maxY; y+=step) {
    cnt = 0;
    for(double x = minX; x <= maxX; x+=step) {
      spc.p[1] = Point(x,y,mcp["xyz"]["z"]);
      spc.trial = spc.p;
      double energy = Energy::systemEnergy(spc,pot,spc.p);
      cnt++;
      
      cout << " " << energy;
      
      // inf 10 5 2
      /*
      if(energy < 1) {
	cout << " " << 0; // No interaction!
      } else if(energy < 3) {
	cout << " " << 1; // capsphereAndSphereCollide No collisoion!
      } else if(energy < 6) {
	cout << " " << 2; // capsphereAndSphereCollide Collisoion!
      } else if(energy < 11) {
	cout << " " << 3; // No caps in play and collision
      } else if(energy < 13) {
	cout << " " << 4; // Nothing at the moment
      } else {
	cout << " " << 5; // Hardsphere-Hardsphere Collision!
      }
      */
      
    }
    cout << endl;
  }
  
  spc.p[0].update();
  spc.trial = spc.p;
  //FormatPQR::save("confout.pqr", spc.p, spc.geo.len); // final PQR snapshot for VMD etc.
  cout << " " << spc.p[0].x() << " " << spc.p[0].y() << " " << spc.p[0].z() << " " << spc.p[0].cap_center << " " << spc.p[0].cap_radius << " " << spc.p[0].radius << " " << spc.p[0].angle_p << " " << spc.p[0].angle_c;
  cout << " " << spc.p[1].x() << " " << spc.p[1].y() << " " << spc.p[1].z() << " " << spc.p[1].cap_center << " " << spc.p[1].cap_radius << " " << spc.p[1].radius << " " << spc.p[1].angle_p << " " << spc.p[1].angle_c;
  cout << " " << mcp["xyz"]["minX"] << " " << mcp["xyz"]["maxX"] << " " << mcp["xyz"]["minY"] << " " << mcp["xyz"]["maxY"] << " " << mcp["xyz"]["step"];
  for(int i = 0; i < cnt-21; i++) {
    cout << " " << 0;
  }
  
  return 0;
  
  
  
  
  
  /*
  spc.p[0] = Point(-5.00001,0,0);
  spc.p[0].cap_center_point = Point(spc.p[0].radius,0,0);
  spc.p[1] = Point(14.9999,0,0);
  spc.p[1].cap_center_point = Point(-spc.p[0].radius,0,0);
  spc.p[2] = Point(0.0,0,0);
  spc.trial = spc.p;
  */
  
  
  

  /*
  spc.p[0] = Point(0,0,0);
  spc.p[0].cap_center_point = Point(spc.p[0].radius,0,0);
  spc.p[0].cap_center = spc.p[0].radius;
  spc.p[0].cap_radius = spc.p[0].radius;
  spc.p[0].is_sphere = false;
  spc.p[1] = Point(15,0,0);
  spc.p[1].cap_center_point = Point(0,0,0);
  spc.p[1].cap_center = 0.0;
  spc.p[1].cap_radius = 0.0;
  spc.p[1].is_sphere = true;
  spc.trial = spc.p;
  
  
  spc.p[1] = Point(5.00001,0,0);
  spc.trial = spc.p;
  for(unsigned int i = 0; i < spc.p.size(); i++)
    cout << "Particle " << i << ": " << atom[spc.p[i].id].name << ", " << spc.p[i].transpose() << endl;
  cout << "Energy: " << Energy::systemEnergy(spc,pot,spc.p) << endl;
  
  spc.p[1] = Point(4.99999,0,0);
  spc.trial = spc.p;
  for(unsigned int i = 0; i < spc.p.size(); i++)
    cout << "Particle " << i << ": " << atom[spc.p[i].id].name << ", " << spc.p[i].transpose() << endl;
  cout << "Energy: " << Energy::systemEnergy(spc,pot,spc.p) << endl;
  
  return 0;
  */
  
  
  spc.load("state");                               // load old config. from disk (if any)
  
  
  sys.init( Energy::systemEnergy(spc,pot,spc.p)  );// store initial total system energy
  
  cout << atom.info() + spc.info() + pot.info() + textio::header("MC Simulation Begins!");
  
  while ( loop[0] ) {  // Markov chain 
    while ( loop[1] ) {
      sys += mv.move();
      analyzer.sample();
      capa.sample(spc);
      xtc.setbox( spc.geo.len );
      xtc.save("traj.xtc", spc.p);
    } // end of micro loop

    sys.checkDrift(Energy::systemEnergy(spc,pot,spc.p)); // compare energy sum with current
    cout << loop.timing();

  } // end of macro loop

  // save to disk
  FormatPQR::save("confout.pqr", spc.p, spc.geo.len); // final PQR snapshot for VMD etc.
  rdf.save("rdf.dat");                // g(r) - not normalized!
  spc.save("state");                     // final simulation state
  capa.save();

  // perform unit tests (irrelevant for the simulation)
  UnitTest test(mcp);                    // class for unit testing
  mv.test(test);
  sys.test(test);

  // print information
  cout << loop.info() + sys.info() + mv.info() + analyzer.info() + test.info() + pot.info();

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
