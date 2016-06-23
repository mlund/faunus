#include <faunus/faunus.h>
#include <faunus/ewald.h>

using namespace Faunus;
using namespace Faunus::Potential;

#define EWALD

typedef Space<Geometry::Cuboid> Tspace;
#ifdef EWALD
typedef LennardJonesLB Tpairpot;
#else
typedef CombinedPairPotential<CoulombWolf,LennardJonesLB> Tpairpot;
#endif

int main() {
  cout << textio::splash();      // show faunus banner and credits
  InputMap mcp("water2.json");   // read input file
  Tspace spc(mcp);
  spc.load("state"); // load old config. from disk (if any)

  // Energy functions and space
#ifdef EWALD
  auto pot = Energy::NonbondedEwald<Tspace,Tpairpot>(mcp)
    + Energy::ExternalPressure<Tspace>(mcp);
#else
  auto pot = Energy::Nonbonded<Tspace,Tpairpot>(mcp)
    + Energy::ExternalPressure<Tspace>(mcp);
#endif
    
  auto pot0 = Energy::NonbondedVector<Tspace,LennardJonesLB>(mcp);
  auto potE = Energy::ExternalPressure<Tspace>(mcp);
  
  auto waters = spc.findMolecules("water");

  // Markov moves and analysis
  Move::Propagator<Tspace> mv( mcp, pot, spc );
  Analysis::RadialDistribution<> rdfOO(0.05);
  Analysis::RadialDistribution<> rdfOH(0.05);
  Analysis::RadialDistribution<> rdfHH(0.05);
  Analysis::MultipoleAnalysis dian(spc,mcp);
  Table2D<double,Average<double> > mucorr;
  mucorr.setResolution(0.05);
  FormatXTC xtc(1000);
  vector<double> energy_T;
  vector<double> energy_LJ;
  vector<double> energy_E;

  EnergyDrift sys;   // class for tracking system energy drifts
#ifdef EWALD
  pot.setSpace(spc);
#endif
  sys.init( Energy::systemEnergy(spc,pot,spc.p)  ); // store total energy

  cout << atom.info() + spc.info() + textio::header("MC Simulation Begins!");
  
  MCLoop loop(mcp);    // class for handling mc loops
  while ( loop[0] ) {          // Markov chain 
    while ( loop[1] ) {

      sys+=mv.move();

      double rnd = slump();
      if ( rnd>0.9 ) {
        rdfOO.sample( spc, atom["OW"].id, atom["OW"].id ); // O-O g(r)
        rdfOH.sample( spc, atom["OW"].id, atom["HW"].id ); // O-H g(r)
        rdfHH.sample( spc, atom["HW"].id, atom["HW"].id ); // H-H g(r)
        dian.sample(spc);
        if ( rnd > 0.99 ) {
          xtc.setbox( spc.geo.len );
          xtc.save( "traj.xtc", spc.p );
        }
	auto beg=spc.groupList().begin();
	auto end=spc.groupList().end();
	for (auto gi=beg; gi!=end; ++gi)
	for (auto gj=gi; ++gj!=end;) {
	  Point mu_a = dipoleMoment(spc,*(*gi));
	  Point mu_b = dipoleMoment(spc,*(*gj));
	  double sca = mu_a.dot(mu_b);
	  double r = spc.geo.dist((*gi)->cm,(*gj)->cm); 
	  mucorr(r) += sca;
	}
	
        energy_T.push_back(Energy::systemEnergy(spc,pot,spc.p));
	energy_LJ.push_back(Energy::systemEnergy(spc,pot0,spc.p));
	energy_E.push_back(Energy::systemEnergy(spc,potE,spc.p));
      }
    } // end of micro loop

    sys.checkDrift(Energy::systemEnergy(spc,pot,spc.p)); // energy drift?

    cout << loop.timing();
  } // end of macro loop

  // save to disk
  FormatPQR::save("confout.pqr", spc.p, spc.geo.len);
  spc.save("state");
  rdfOO.save("rdfOO.dat");
  rdfOH.save("rdfOH.dat");
  rdfHH.save("rdfHH.dat");
  mucorr.save("mucorr.dat");
  dian.save();
  spc.save("state");
  
  string file = "energy.dat";
  std::ofstream f(file.c_str());
  if (f)
    for (unsigned int i = 0; i < energy_T.size(); i++)
      f << std::left << std::setw(10) << energy_T.at(i) << "     " << energy_LJ.at(i) << "      " << energy_E.at(i) << endl;

  // perform unit 
  UnitTest test(mcp);
  sys.test(test);
  mv.test(test);

  // print information
  cout << loop.info() + pot.info() + sys.info() + mv.info() + dian.info() + test.info();

  return test.numFailed();
}

/**  @page example_water2 Example: SPC Water (V2)

 This will simulate SPC water in a cubic volume using
 the Wolf method for electrostatic interactions.
 This version uses a fake cell list to discard
 interactions beyond a specified water-water mass-center
 cutoff.

 Run this example from the main faunus directory:

 ~~~~~~~~~~~~~~~~~~~
 $ make example_water2
 $ cd src/examples
 $ ./water2.run
 ~~~~~~~~~~~~~~~~~~~

 ![Water](water.png)

 water2.cpp
 ==========

 @includelineno examples/water2.cpp

 water2.json
 -----------

 @includelineno examples/water2.json

*/
