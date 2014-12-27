#include <faunus/faunus.h>
using namespace Faunus;

#ifdef CUBOID
typedef Geometry::Cuboid Tgeometry;   // specify geometry - here cube w. periodic boundaries
#else
typedef Geometry::Sphere Tgeometry;   // sphere with hard boundaries
#endif
typedef Potential::CoulombHS Tpairpot;// particle pair potential: primitive model

typedef Space<Tgeometry,PointParticle> Tspace;

int main() {
  cout << textio::splash();           // show faunus banner and credits
  
  InputMap mcp("polymers.input");     // open user input file
  MCLoop loop(mcp);                   // class for handling mc loops
  EnergyDrift sys;                    // class for tracking system energy drifts
  UnitTest test(mcp);                 // class for unit testing

  Tspace spc(mcp);                    // instantiate space and load molecules/atoms

  // Create Space and a Hamiltonian with three terms
  auto pot = Energy::Nonbonded<Tspace,Tpairpot>(mcp)
    + Energy::ExternalPressure<Tspace>(mcp) + Energy::Bonded<Tspace>();

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

  // locate salt and polymers
  auto pol = spc.findMolecules("polymer");               // grop vector w. all molecules named "polymer" 
  Group all = Group(0, spc.p.size()-1);                  // group w. all atoms

  cout << atom.info() << spc.info() << pot.info() << textio::header("MC Simulation Begins!");

  while ( loop[0] ) {  // Markov chain 
    while ( loop[1] ) {
      int k, i = slump.range(0,5);
      switch (i) {
        case 0:
          mv.setGroup(all);
          sys+=mv.move( all.size() );    // translate salt and monomers
          for (auto &g : pol) {
            g->setMassCenter(spc);
            shape.sample(*g,spc);        // sample gyration radii etc.
          }
          break;
        case 2:
          k=pol.size();
          while (k-->0) {
            gmv.setGroup( **slump.element(pol.begin(), pol.end() ) );
            sys+=gmv.move();            // translate/rotate polymers
          }
          break;
        case 3:
          sys+=iso.move();              // isobaric volume move
          break;
        case 4:
          k=pol.size();
          while (k-->0) {
            crank.setGroup( **slump.element(pol.begin(), pol.end() ) );
            sys+=crank.move();          // crank shaft move
          }
          break;
        case 5:
          k=pol.size();
          while (k-->0) {
            pivot.setGroup( **slump.element(pol.begin(), pol.end() ) );
            sys+=pivot.move();          // pivot move
          }
          break;
      }

      // polymer-polymer mass center rdf
      for (auto i = pol.begin(); i != pol.end(); ++i )
        for (auto j = i; ++j != pol.end(); )
          rdf( spc.geo.dist( (**i).cm, (**j).cm) )++;

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
 * $ ./polymers.run
 * ~~~~~~~~~~~~~~~~~~~
 *
 * polymers.cpp
 * ============
 *
 * @includelineno examples/polymers.cpp
 */

