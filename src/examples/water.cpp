#include <faunus/faunus.h>
using namespace Faunus;
using namespace Faunus::Potential;

typedef Geometry::Cuboid Tgeometry; // specify geometry 
typedef CombinedPairPotential<CoulombWolf,LennardJonesLB> Tpairpot;
typedef Space<Tgeometry,PointParticle> Tspace;

int main() {
  cout << textio::splash();         // show faunus banner and credits
                                  
  InputMap mcp("water.input");      // open user input file
  MCLoop loop(mcp);                 // class for handling mc loops
  EnergyDrift sys;                  // class for tracking system energy drifts
  UnitTest test(mcp);               // class for unit testing

  // Create Space and a Hamiltonian (nonbonded+NVT)
  Tspace spc(mcp);
  auto pot = Energy::Nonbonded<Tspace,Tpairpot>(mcp)
    + Energy::ExternalPressure<Tspace>(mcp);

  // Read single water from disk and add N times
  Group sol;
  sol.setMolSize(3);
  string file = mcp.get<string>("mol_file");
  int N=mcp("mol_N",1);
  for (int i=0; i<N; i++) {
    Tspace::ParticleVector v;
    FormatAAM::load(file,v);
    Geometry::FindSpace().find(spc.geo, spc.p, v);
    Group g = spc.insert(v);// Insert into Space
    sol.setrange(0, g.back());
  }
  spc.enroll(sol);

  // Markov moves and analysis
  Move::Isobaric<Tspace> iso(mcp,pot,spc);
  Move::TranslateRotate<Tspace> gmv(mcp,pot,spc);
  Analysis::RadialDistribution<> rdf(0.05);
  Analysis::DielectricConstant gdc(spc);

  spc.load("state");                               // load old config. from disk (if any)
  sys.init( Energy::systemEnergy(spc,pot,spc.p)  );// store init system energy

  cout << atom.info() + spc.info() + pot.info() + textio::header("MC Simulation Begins!");

  while ( loop.macroCnt() ) {                      // Markov chain 
    while ( loop.microCnt() ) {
      int j,i=slp_global.rand() % 2;
      int k=sol.numMolecules();                    //number of water molecules
      Group g;
      switch (i) {
        case 0:
          while (k-->0) {
            j=sol.randomMol();                     // pick random water mol.
            sol.getMolecule(j,g);
            g.name="water";
            g.setMassCenter(spc);                  // mass center needed for rotation
            gmv.setGroup(g);
            sys+=gmv.move();                       // translate/rotate
          }
          break;
        case 1:
          sys+=iso.move();                         // volume move
          break;
      }
      gdc.samplePP(spc.geo,spc);
      // sample oxygen-oxygen rdf
      if (slp_global()>0.9) {
        auto id = atom["OW"].id;
        rdf.sample(spc,sol,id,id);
      }
    
    } // end of micro loop
    //cout << gdc.info() << endl;
    //cout << "Inf: " << gdc.getDielInfty() << endl;
    sys.checkDrift(Energy::systemEnergy(spc,pot,spc.p));
    cout << loop.timing();

  } // end of macro loop

  //cout << gdc.info() << endl;
  rdf.save("rdf.dat");
  spc.save("state");
  FormatPQR::save("confout.pqr", spc.p);

  // perform unit tests
  iso.test(test);
  gmv.test(test);
  sys.test(test);

  // print information
  cout << loop.info() + sys.info() + gmv.info() + iso.info() + test.info();

  return test.numFailed();
}
/**  @page example_water Example: SPC Water
 *
 This will simulate SPC water in a cubic volume using
 the Wolf method for electrostatic interactions.

 Run this example from the `examples` directory:

 ~~~~~~~~~~~~~~~~~~~
 $ make
 $ cd src/examples
 $ ./water.run
 ~~~~~~~~~~~~~~~~~~~

 water.cpp
 ============

 @includelineno examples/water.cpp

*/

