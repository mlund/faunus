#include <faunus/faunus.h>
using namespace Faunus;

#include "power_sasa.h"
#include <time.h>


#ifdef CUBOID
typedef Geometry::Cuboid Tgeometry;   // specify geometry - here cube w. periodic boundaries
#else
typedef Geometry::Sphere Tgeometry;   // sphere with hard boundaries
#endif
typedef Potential::HardSphere Tpairpot;// particle pair potential: primitive model

typedef Space<Tgeometry,PointParticle> Tspace;

//typedef Eigen::Vector3f Coord;
typedef Point Coord;
typedef float Scalar;

int main() {
  
  std::vector<Coord> sasaCoords;     //vector of sasa coordinates
  std::vector<Scalar> sasaWeights;   //vector of sasa weights
  
  clock_t clk_a, clk_b, clk_sum = 0;
  
  /*
  sasaCoords.push_back(Coord(1.0,2.0,3.0));
  sasaWeights.push_back(1.0);
  
  POWERSASA::PowerSasa<Scalar,Coord> *ps = 
    new POWERSASA::PowerSasa<Scalar,Coord>(sasaCoords, sasaWeights, 1, 1, 1, 1);
  ps->calc_sasa_all();
  
  
  Scalar volume = 0.0, surf = 0.0;
  for (unsigned int i = 0; i < sasaCoords.size(); ++i)
  {
	printf("%4d sasa=%7.3lf vol=%7.3lf\n", i, (ps->getSasa())[i], (ps->getVol())[i]);
	volume += (ps->getVol())[i];
	surf += (ps->getSasa())[i];
  }
  printf("volume=%lf\n", volume);
  printf("sasa  =%lf\n", surf);

  delete ps;
  */
  cout << textio::splash();           // show faunus banner and credits
  
  InputMap mcp("polymers.input");     // open user input file
  MCLoop loop(mcp);                   // class for handling mc loops
  EnergyDrift sys;                    // class for tracking system energy drifts
  UnitTest test(mcp);                 // class for unit testing

  // Create Space and a Hamiltonian with three terms
  Tspace spc(mcp);
  auto pot = Energy::Nonbonded<Tspace,Tpairpot>(mcp)
    + Energy::SASAEnergy<Tspace>(mcp)
    + Energy::ExternalPressure<Tspace>(mcp) + Energy::Bonded<Tspace>();

  auto bonded = &pot.second; // pointer to bond energy class

  // Add salt
  Group salt;
  salt.addParticles(spc, mcp);

  // Add polymers
  vector<Group> pol( mcp.get("polymer_N",0));            // vector of polymers
  string file = mcp.get<string>("polymer_file", "");
  double req = mcp.get<double>("polymer_eqdist", 0);
  double k   = mcp.get<double>("polymer_forceconst", 0);
  for (auto &g : pol) {                                  // load polymers
    Tspace::ParticleVector v;                            // temporary, empty particle vector
    FormatAAM::load(file,v);                             // load AAM structure into v
    Geometry::FindSpace().find(spc.geo,spc.p,v);         // find empty spot in particle vector
    g = spc.insert(v);                                   // insert into space
    g.name="Polymer";
    spc.enroll(g);
    for (int i=g.front(); i<g.back(); i++)
      bonded->add(i, i+1, Potential::Harmonic(k,req));   // add bonds
  }
  Group allpol( pol.front().front(), pol.back().back() );// make group w. all polymers

  // Markov moves and analysis
  Move::Isobaric<Tspace> iso(mcp,pot,spc);
  Move::TranslateRotate<Tspace> gmv(mcp,pot,spc);
  Move::AtomicTranslation<Tspace> mv(mcp, pot, spc);
  Move::CrankShaft<Tspace> crank(mcp, pot, spc);
  Move::Pivot<Tspace> pivot(mcp, pot, spc);
  Analysis::PolymerShape shape;
  Analysis::RadialDistribution<> rdf(0.2);
  Scatter::DebyeFormula<Scatter::FormFactorUnity<>> debye(mcp);

  spc.load("state");                                     // load old config. from disk (if any)
  sys.init( Energy::systemEnergy(spc,pot,spc.p)  );      // store initial total system energy

  cout << atom.info() << spc.info() << pot.info() << textio::header("MC Simulation Begins!");

  while ( loop.macroCnt() ) {  // Markov chain 
    while ( loop.microCnt() ) {
      int k,i=slp_global.rand() % 6;
      switch (i) {
        case 0:
          mv.setGroup(salt);
          sys+=mv.move( salt.size() );  // translate salt
          break;
        case 1:
          mv.setGroup(allpol);
          sys+=mv.move( allpol.size() );// translate monomers
          for (auto &g : pol) {
            g.setMassCenter(spc);
            shape.sample(g,spc);        // sample gyration radii etc.
          }
          break;
        case 2:
          k=pol.size();
          while (k-->0) {
            gmv.setGroup( pol[ slp_global.rand() % pol.size() ] );
            sys+=gmv.move();            // translate/rotate polymers
          }
          break;
        case 3:
          sys+=iso.move();              // isobaric volume move
          break;
        case 4:
          k=pol.size();
          while (k-->0) {
            crank.setGroup( pol[ slp_global.rand() % pol.size() ] );
            sys+=crank.move();          // crank shaft move
          }
          break;
        case 5:
          k=pol.size();
          while (k-->0) {
            pivot.setGroup( pol[ slp_global.rand() % pol.size() ] );
            sys+=pivot.move();          // pivot move
          }
          break;
      }

      // polymer-polymer mass center rdf
      for (auto i=pol.begin(); i!=pol.end()-1; i++)
        for (auto j=i+1; j!=pol.end(); j++)
          rdf( spc.geo.dist(i->cm,j->cm) )++;


      // SASA
      sasaCoords.clear();
      sasaWeights.clear();

/*
      for (auto i = pol.begin(); i != pol.end(); i++) {
		sasaCoords.push_back(i->cm);
		sasaWeights.push_back(5.0);
	  }
*/

      /*
      for (auto i = allpol.begin(); i != allpol.end(); i++) {
		auto monomer = spc.p.at(*i);
		sasaCoords.push_back(monomer);
		sasaWeights.push_back(10.0);
	  }
      clk_a = clock();
      auto ps = new POWERSASA::PowerSasa<Scalar,Coord>(sasaCoords, sasaWeights, 1, 1, 1, 1);
      ps->calc_sasa_all();
      clk_b = clock();
      clk_sum += clk_b - clk_a;
      
      {
		  Scalar volume = 0.0, surf = 0.0;
		  for (unsigned int i = 0; i < sasaCoords.size(); ++i)
		  {
			volume += (ps->getVol())[i];
			surf += (ps->getSasa())[i];
		  }
		  //printf("volume=%lf\n", volume);
		  //printf("sasa  =%lf\n", surf);
	  }
      delete ps;
*/
    } // end of micro loop

    // sample scattering
    debye.sample(spc.p);

    sys.checkDrift(Energy::systemEnergy(spc,pot,spc.p)); // compare energy sum with current
    cout << loop.timing();

  } // end of macro loop

  // save to disk
  rdf.save("rdf_p2p.dat");
  spc.save("state");
  debye.save("debye.dat");
  FormatPQR::save("confout.pqr", spc.p);

  // perform unit tests
  iso.test(test);
  gmv.test(test);
  mv.test(test);
  sys.test(test);

  // print information
  cout << loop.info() << sys.info() << mv.info() << gmv.info() << crank.info()
    << pivot.info() << iso.info() << shape.info() << test.info();

  cout << endl << endl << "SASA" << endl;
  cout << "Time " << (float) clk_sum/CLOCKS_PER_SEC << " s = " << (float) clk_sum/CLOCKS_PER_SEC/3600 << " h" << endl;
  
  cout << sasaCoords.size() << endl;

  cout << pot.first.first.second.info() ;
  
  /*
  cout << allpol.size() << endl;
  cout << allpol.info() << endl;
  cout << spc.info() << endl;
  cout << typeid(spc.p.at(34)).name() << endl;
  cout << typeid(allpol.begin()).name() << endl;
  cout << spc.findGroup(34)->info() << endl;
*/	
  
  //for (auto i = allpol.begin(); i != allpol.end(); i++) {
  //}

  return test.numFailed();
}
/*! \page example_polymers Example: Polymers
 *
 This will simulate an arbitrary number of linear polymers in the NPT ensemble
 with explicit salt particles and implicit solvent (dielectric continuum).
 We include the following Monte Carlo moves:

 - salt translation
 - monomer translation
 - polymer translation and rotation
 - polymer crankshaft and pivot rotations
 - isobaric volume move (NPT ensemble)

 ![Hardsphere polyelectrolytes with counter ions](polymers.png)

 Run this example from the `examples` directory:

 ~~~~~~~~~~~~~~~~~~~
 $ make
 $ cd src/examples
 $ ./polymers.run
 ~~~~~~~~~~~~~~~~~~~

 polymers.cpp
 ============

 \includelineno examples/polymers.cpp

*/

