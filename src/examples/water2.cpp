#include <faunus/faunus.h>

using namespace Faunus;
using namespace Faunus::Potential;

typedef Space<Geometry::Cuboid> Tspace;
typedef CombinedPairPotential<CoulombWolf,LennardJonesLB> Tpairpot;

int main() {

  cout << textio::splash();      // show faunus banner and credits
  InputMap mcp("water2.input");  // read input file

  // Energy functions and space
  auto pot = Energy::NonbondedCutg2g<Tspace,Tpairpot>(mcp)
    + Energy::ExternalPressure<Tspace>(mcp);
  Tspace spc(mcp);

  // Load and add molecules to Space
  auto N    = mcp.get<int>("mol_N",1);
  auto file = mcp.get<string>("mol_file");
  vector<Group> water(N);
  Tspace::ParticleVector v;      // temporary, empty particle vector
  FormatAAM::load(file,v);       // load AAM structure into v
  for (auto &i : water) {
    Geometry::FindSpace f;
    f.allowMatterOverlap=true;
    f.find(spc.geo,spc.p,v);     // find empty spot in particle vector
    i = spc.insert(v);           // Insert into Space
    i.name="h2o";
    spc.enroll(i);
  }

  // Markov moves and analysis
  Move::TranslateRotate<Tspace> gmv(mcp,pot,spc);
  Move::Isobaric<Tspace> iso(mcp,pot,spc);
  Analysis::RadialDistribution<> rdf(0.05);

  spc.load("state"); // load old config. from disk (if any)

  EnergyDrift sys;   // class for tracking system energy drifts
  sys.init( Energy::systemEnergy(spc,pot,spc.p)  ); // store total energy

  cout << atom.info() + spc.info() + pot.info() + textio::header("MC Simulation Begins!");

  MCLoop loop(mcp);            // class for handling mc loops
  while ( loop.macroCnt() ) {  // Markov chain 
    while ( loop.microCnt() ) {
      int i=slp_global.rand() % 2;
      int j,k=water.size();
      Group g;
      switch (i) {
        case 0:
          while (k-->0) {
            j=slp_global.rand() % (water.size());
            gmv.setGroup(water[j]);
            sys+=gmv.move();   // translate/rotate polymers
          }
          break;
        case 1:
          sys+=iso.move();
          break;
      }

      // sample oxygen-oxygen rdf
      if (slp_global()>0.9) {
        auto id = atom["OW"].id;
        rdf.sample(spc,id,id);
      }

    } // end of micro loop

    sys.checkDrift(Energy::systemEnergy(spc,pot,spc.p)); // energy drift?

    cout << loop.timing();
  } // end of macro loop

  // save to disk
  FormatPQR::save("confout.pqr", spc.p);
  spc.save("state");
  rdf.save("rdf.dat");
  spc.save("state");
  FormatPQR::save("confout.pqr", spc.p, spc.geo.len);

  // perform unit 
  UnitTest test(mcp);
  iso.test(test);
  gmv.test(test);
  sys.test(test);

  // print information
  cout << loop.info() + sys.info() + gmv.info() + iso.info() + test.info();
}

/**  @page example_water2 Example: SPC Water (V2)
 *
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
 ============

 @includelineno examples/water2.cpp

*/
