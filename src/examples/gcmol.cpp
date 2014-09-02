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
  
  InputMap mcp("gcmol.input");     // open user input file
  MCLoop loop(mcp);                   // class for handling mc loops
  EnergyDrift sys;                    // class for tracking system energy drifts
  UnitTest test(mcp);

  // Create Space and a Hamiltonian with three terms
  Tspace spc(mcp);
  auto pot = Energy::Nonbonded<Tspace,Tpairpot>(mcp);

  string file = mcp.get<string>("polymer_file", "");


  //
  //    ADDING SPECIES, NOTE: MUST ALLOCATE ON HEAP
  //
/*"topology": {
    "salt": { "atoms": "Na,Cl", "init":"RANDOM"},
    "salt2": { "atoms": "Mg,Cl,Cl", "init":"RANDOM"},
    "chloride": { "atoms": "Cl,Cl,Cl,Cl", "init":"RANDOM"},
    "polymer": {"activity": 0.000000000000000005, "atoms": "MM,MM,MM,MM", "init": "POOL"},
    "polymer2": {"activity": 0.000000000000000005, "atoms": "MM,MM,MM,MM", "init": "POOL"}
  }*/

  string na("Na");
  string cl("Cl");
  string mg("Mg");

  Group *salt = new Group("salt"); // group for salt particles
  salt->setfront(spc.p.size());
  salt->addParticles(spc, na, 8);
  salt->addParticles(spc, cl, 8);
  salt->setMassCenter(spc);
  spc.enroll(*salt);

  // Add polymers
  vector<Group*> pol( mcp.get("polymer_N",0));            // vector of polymers

  for (auto* g : pol) { // load polymers
      g = new Group();
    Tspace::ParticleVector v;                            // temporary, empty particle vector
    FormatAAM::load(file,v);                             // load AAM structure into v
    Geometry::FindSpace().find(spc.geo,spc.p,v);         // find empty spot in particle vector
    *g = spc.insert(v);
    g->name="polymer";
    spc.enroll(*g);
  }

  Group *salt2 = new Group("salt2"); // group for salt particles
  salt2->setfront(spc.p.size());
  salt2->addParticles(spc, mg, 8);
  salt2->addParticles(spc, cl, 16);
  salt2->setMassCenter(spc);
  spc.enroll(*salt2);

  // Add polymers
  vector<Group*> pol2( mcp.get("polymer2_N",0));            // vector of polymers

  for (auto* g : pol2) { // load polymers
      g = new Group();
    Tspace::ParticleVector v;                            // temporary, empty particle vector
    FormatAAM::load(file,v);                             // load AAM structure into v
    Geometry::FindSpace().find(spc.geo,spc.p,v);         // find empty spot in particle vector
    *g = spc.insert(v);
    g->name="polymer2";
    spc.enroll(*g);
  }

  Group *chloride = new Group("chloride"); // group for salt particles
  chloride->setfront(spc.p.size());
  chloride->addParticles(spc, cl, 16);
  chloride->setMassCenter(spc);
  spc.enroll(*chloride);

  salt = salt2 = chloride = NULL;

  //
  //    ALL SPECIES ADDED
  //

  // Markov moves and analysis
  Move::GCMolecular<Tspace> gc(mcp, pot, spc);
  Move::TranslateRotate<Tspace> mv(mcp,pot,spc);

  //////////////////////////////////////////////////////////////////////////////////////
  //
  //    Configuration = p_vec is one molecule large
  //
  p_vec conf;
  FormatAAM::load(file,conf);
  molecule.pushConfiguration("polymer", conf);
  molecule.pushConfiguration("polymer2",conf);
  //////////////////////////////////////////////////////////////////////////////////////

  //spc.load("state");                                     // load old config. from disk (if any)
  sys.init( Energy::systemEnergy(spc,pot,spc.p)  );      // store initial total system energy

  cout << atom.info() << molecule.info() << gc.info() << spc.info() << pot.info() << textio::header("MC Simulation Begins!");

    while ( loop.macroCnt() ) {  // Markov chain
        while ( loop.microCnt() ) {
            sys+=gc.move();
        }
        sys.checkDrift(Energy::systemEnergy(spc,pot,spc.p)); // compare energy sum with current
        cout << loop.timing();
    } // end of macro loop

  // save to disk

  spc.save("state");
  FormatPQR::save("confout.pqr", spc.p);

  sys.test(test);
  gc.test(test);

  // print information
  cout << loop.info() << sys.info() << gc.info() << spc.info() << test.info();

  //
  // clean allocated memory
  //
  spc.freeGroups();

  return test.numFailed();
}
