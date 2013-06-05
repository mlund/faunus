#include <faunus/faunus.h>
#include <faunus/multipole.h>
using namespace Faunus;
using namespace Faunus::Move;
using namespace Faunus::Potential;

typedef Space<Geometry::Cuboid,DipoleParticle> Tspace;
typedef CombinedPairPotential<LennardJones,DipoleDipole> Tpairpot;
//typedef CombinedPairPotential<Tpairpot1,IonDipole> Tpairpot2;
//typedef CombinedPairPotential<Tpairpot2,IonQuad> Tpairpot;

int main() {
  cout << textio::splash();         // show faunus banner and credits
                                  
  InputMap mcp("water.input");      // open user input file
  MCLoop loop(mcp);                 // class for handling mc loops
  EnergyDrift sys;                  // class for tracking system energy drifts
  UnitTest test(mcp);               // class for unit testing

  // Create Space and a Hamiltonian (nonbonded+NVT)
  Tspace spc(mcp);
  auto pot = Energy::NonbondedVector<Tspace,Tpairpot>(mcp)
    + Energy::ExternalPressure<Tspace>(mcp) + Energy::Bonded<Tspace>();

    
      auto bonded = &pot.second; // pointer to bond energy class
  // Read single water from disk and add N times
  Group sol;
  sol.setMolSize(5);
  string file = mcp.get<string>("mol_file");
  int N=mcp("mol_N",1);
  for (int i=0; i<N; i++) {
    Tspace::ParticleVector v;
    FormatAAM::load(file,v);
    Geometry::FindSpace().find(spc.geo, spc.p, v);
    Group g = spc.insert(v);// Insert into Space
    sol.setrange(0, g.back());
    bonded->add(g.front(), g.front()+3, Potential::Harmonic(0.557,0.4555));       // add bonds O -B1
    bonded->add(g.front(), g.front()+4, Potential::Harmonic(0.557,0.4555));       // add bonds O -B2
    bonded->add(g.front()+1, g.front()+3, Potential::Harmonic(0.557,0.4555));     // add bonds H1-B1
    bonded->add(g.front()+2, g.front()+4, Potential::Harmonic(0.557,0.4555));     // add bonds H2-B2
    bonded->add(g.front()+1, g.front()+2, Potential::Harmonic(0.557,1.008568));   // add bonds H1-H2
  }
  spc.enroll(sol);

  // Markov moves and analysis
  Move::Isobaric<Tspace> iso(mcp,pot,spc);
  Move::TranslateRotate<Tspace> gmv(mcp,pot,spc);
  Analysis::RadialDistribution<> rdf(0.05);
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

      // sample oxygen-oxygen rdf
      if (slp_global()>0.9) {
        auto id = atom["OW"].id;
        rdf.sample(spc,sol,id,id);
      }

    } // end of micro loop

    sys.checkDrift(Energy::systemEnergy(spc,pot,spc.p));
    cout << loop.timing();

  } // end of macro loop

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
 This will simulate an arbitrary SPC water in a cubic box using
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

