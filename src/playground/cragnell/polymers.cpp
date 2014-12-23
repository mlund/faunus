#include <faunus/faunus.h>
using namespace Faunus;

typedef Geometry::Cuboid Tgeometry;            // specify geometry - here cube w. periodic boundaries
typedef Potential::DebyeHuckelLJ Tpairpot;     // particle pair potential: primitive model
typedef Space<Tgeometry,PointParticle> Tspace;

int main() {
  cout << textio::splash();           // show faunus banner and credits
  
  InputMap mcp("polymers.input");     // open user input file
  MCLoop loop(mcp);                   // class for handling mc loops
  EnergyDrift sys;                    // class for tracking system energy drifts

  // Create Space and a Hamiltonian with three terms
  Tspace spc(mcp);
  auto pot = Energy::Nonbonded<Tspace,Tpairpot>(mcp)
    + Energy::Bonded<Tspace>(); //+ Energy::EquilibriumEnergy<Tspace>(mcp);

  auto bonded = &pot.second; // pointer to bond energy class

  // Add salt
  Group salt;
  salt.addParticles(spc, mcp);

  // Add polymers
  vector<Group> pol( mcp.get("polymer_N",0));            // vector of polymers
  string file = mcp.get<string>("polymer_file", "");
  double req = mcp.get<double>("polymer_eqdist", 0);
  double k   = mcp.get<double>("polymer_forceconst", 0);

  auto bondpot = Potential::Harmonic(k,req) - Potential::LennardJones(mcp, "blj");

  for (auto &g : pol) {                                  // load polymers
    Tspace::ParticleVector v;                            // temporary, empty particle vector
    FormatAAM::load(file,v);                             // load AAM structure into v
    Geometry::FindSpace().find(spc.geo,spc.p,v);         // find empty spot in particle vector
    g = spc.insert(v);                                   // insert into space
    g.name="Polymer";
    spc.enroll(g);
    for (int i=g.front(); i<g.back(); i++)
      bonded->add(i, i+1, bondpot);   // add bonds
  }
  Group allpol( pol.front().front(), pol.back().back() );// make group w. all polymers

  // Markov moves and analysis
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
      int k,i= slump.rand() % 5;
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
            gmv.setGroup( pol[ slump.rand() % pol.size() ] );
            sys+=gmv.move();            // translate/rotate polymers
          }
          break;
        case 3:
          k=pol.size();
          while (k-->0) {
            crank.setGroup( pol[ slump.rand() % pol.size() ] );
            sys+=crank.move();          // crank shaft move
          }
          break;
        case 4:
          k=pol.size();
          while (k-->0) {
            pivot.setGroup( pol[ slump.rand() % pol.size() ] );
            sys+=pivot.move();          // pivot move
          }
          break;
      }

      // polymer-polymer mass center rdf
      for (auto i=pol.begin(); i!=pol.end()-1; i++)
        for (auto j=i+1; j!=pol.end(); j++)
          rdf( spc.geo.dist(i->cm,j->cm) )++;

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

  // print information
  cout << loop.info() << sys.info() << mv.info() << gmv.info() << crank.info()
    << pivot.info() << shape.info();
}
