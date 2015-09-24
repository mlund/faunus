#include <faunus/faunus.h>
using namespace Faunus;

typedef Geometry::Cuboid Tgeometry;   // specify geometry - here cube w. periodic boundaries
typedef Potential::CoulombHS Tpairpot;// particle pair potential: primitive model

typedef Space<Tgeometry,PointParticle> Tspace;

int main() {
  cout << textio::splash();           // show faunus banner and credits
  
  InputMap mcp("polyelectrolytes.json");      // open user input file
  MCLoop loop(mcp);                   // class for handling mc loops
  EnergyDrift sys;                    // class for tracking system energy drifts
  UnitTest test(mcp);                 // class for unit testing

  Tspace spc(mcp);                    // instantiate space and load molecules/atoms
  spc.load("state", Tspace::RESIZE);  // load old config. from disk (if any)

  // Create Space and a Hamiltonian with three terms
  auto pot = Energy::Nonbonded<Tspace,Tpairpot>(mcp)
    + Energy::ExternalPressure<Tspace>(mcp) + Energy::Bonded<Tspace>()
    + Energy::EquilibriumEnergy<Tspace>(mcp);

  Move::Propagator<Tspace> mv( mcp, pot, spc );
  Analysis::CombinedAnalysis analyze(mcp, pot, spc);

  //Analysis::PolymerShape shape;
  Analysis::RadialDistribution<> rdf(0.2);
  //Analysis::ChargeMultipole cm;

  //spc.load("state");                                // load old config. from disk (if any)

  auto pol = spc.findMolecules("polymer");          // all molecules named "polymer" 
  auto salt = spc.findMolecules("salt");
  std::cout << salt.size()<<std::endl;
  if (salt.size() != 1 )
    std::cout << "Number of salt groups different from ONE!!!"<<std::endl;
  pot.first.first.second.ignore.insert(salt[0]);
//  pot.first.second.ignore.insert(salt[0]);
  sys.init( Energy::systemEnergy(spc,pot,spc.p)  ); // store initial total system energy
  cout << atom.info() << spc.info() << pot.info() << textio::header("MC Simulation Begins!");

  while ( loop[0] ) {  // Markov chain 
    while ( loop[1] ) {
      sys += mv.move();

      analyze.sample();
      //for (auto i : pol) {                            // sample gyration radii etc. 
      //  shape.sample(*i,spc);
      //  cm.sample(*i, spc);
      //}

      for (auto i = pol.begin(); i != pol.end(); ++i )// sample cm-cm g(r)
        for (auto j = i; ++j != pol.end(); )
          rdf( spc.geo.dist( (**i).cm, (**j).cm) )++;

    } // end of micro loop

    sys.checkDrift(Energy::systemEnergy(spc,pot,spc.p)); // compare energy sum with current
    cout << loop.timing();

  } // end of macro loop
  std::cout << "Anything?"<<std::endl;
  // save to disk
  rdf.save("rdf_p2p.dat");
  spc.save("state");
  FormatPQR::save("confout.pqr", spc.p);

  // unit tests
  mv.test(test);
  sys.test(test);

  // print information
  cout << loop.info() << mv.info() << analyze.info() //<< cm.info()
    << sys.info() << test.info();

  return test.numFailed();
}
/** 
 * @page example_polymers Example: Polymers
 *
 * This will simulate an arbitrary number of linear polymers in the NPT ensemble
 * with explicit salt particles and implicit solvent (dielectric continuum).
 * We include the following Monte Carlo moves:
 *
 * - salt translation
 * - monomer translation
 * - polymer translation and rotation
 * - polymer crankshaft and pivot rotations
 * - isobaric volume move (NPT ensemble)
 *
 * ![Hardsphere polyelectrolytes with counter ions](polymers.png)
 *
 * Run this example from the `examples` directory:
 *
 * ~~~~~~~~~~~~~~~~~~~
 * $ make
 * $ cd src/examples
 * $ python ./polyelectrolytes.py
 * ~~~~~~~~~~~~~~~~~~~
 *
 * polyelectrolytes.py
 * ===========
 * This is a python run-script that generates all required
 * input as well as run the main program. This is also used
 * for running the unittest when running `make test`. 
 *
 * @include polyelectrolytes.py
 */

