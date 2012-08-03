/*
   This will simulate an arbitrary number of linear polymers in an an NPT simulation
   container with explicit salt particles and implicit solvent (dielectric continuum).
   We include the following Monte Carlo move:

   1) salt translation
   2) polymer translation and rotation
   3) polymer crankshaft rotation
   4) polymer pivot rotation
   5) monomer translation
   6) isobaric volume move (NPT ensemble)
   */

#include <faunus/faunus.h>

using namespace Faunus;

#ifdef CUBOID
typedef Geometry::Cuboid Tgeometry;   // specify geometry - here a cube with periodic boundaries
#else
typedef Geometry::Sphere Tgeometry;   // here a sphere with hard boundaries
#endif
typedef Potential::CoulombHS Tpairpot;// particle pair potential

int main() {
  InputMap mcp("polymers.input");      // open user input file
  MCLoop loop(mcp);                    // class for handling mc loops
  FormatPQR pqr;                       // PQR structure file I/O
  FormatAAM aam;                       // AAM structure file I/O
  EnergyDrift sys;                     // class for tracking system energy drifts
  UnitTest test(mcp);                  // class for unit testing

  // Energy functions and space
  Energy::Hamiltonian pot;
  auto nonbonded = pot.create( Energy::Nonbonded<Tpairpot,Tgeometry>(mcp) );
  auto bonded    = pot.create( Energy::Bonded() );
  Space spc( pot.getGeometry() );

  // Markov moves and analysis
  Move::Isobaric iso(mcp,pot,spc);
  Move::TranslateRotate gmv(mcp,pot,spc);
  Move::AtomicTranslation mv(mcp, pot, spc);
  Move::CrankShaft crank(mcp, pot, spc);
  Move::Pivot pivot(mcp, pot, spc);
  Analysis::PolymerShape shape;
  Analysis::RadialDistribution<float,int> rdf(0.2);

  // Add salt
  GroupAtomic salt(spc, mcp);
  salt.name="Salt";

  // Add polymers
  vector<GroupMolecular> pol( mcp.get("polymer_N",0));    // vector of polymers
  string polyfile = mcp.get<string>("polymer_file", "");
  double req    = mcp.get<double>("polymer_eqdist", 0);
  double k      = mcp.get<double>("polymer_forceconst", 0);
  atom["MM"].dp = 10.;
  for (auto &g : pol) {                                   // load polymers
    aam.load(polyfile);
    Geometry::FindSpace f;
    f.find(*spc.geo, spc.p, aam.particles() );            // find empty spot in particle vector
    g = spc.insert( aam.particles() );                    // insert into space
    g.name="Polymer";
    spc.enroll(g);
    for (int i=g.front(); i<g.back(); i++)
      bonded->add(i, i+1, Potential::Harmonic(k,req));    // add bonds
  }
  Group allpol( pol.front().front(), pol.back().back() ); // make group w. all polymers

  spc.load("state");                                      // load old configuration from disk (if any)

  sys.init( Energy::systemEnergy(spc,pot,spc.p)  );       // store initial total system energy

  cout << atom.info() << spc.info() << pot.info() << textio::header("MC Simulation Begins!");

  while ( loop.macroCnt() ) {  // Markov chain 
    while ( loop.microCnt() ) {
      int k,i=rand() % 6;
      switch (i) {
        case 0:
          mv.setGroup(salt);
          sys+=mv.move( salt.size() ); // translate salt
          break;
        case 1:
          mv.setGroup(allpol);
          sys+=mv.move( allpol.size() ); // translate monomers
          for (auto &g : pol) {
            g.setMassCenter(spc);
            shape.sample(g,spc);
          }
          break;
        case 2:
          k=pol.size();
          while (k-->0) {
            gmv.setGroup( pol[ rand() % pol.size() ] );
            sys+=gmv.move(); // translate/rotate polymers
          }
          break;
        case 3:
          sys+=iso.move(); // isobaric volume move
          break;
        case 4:
          k=pol.size();
          while (k-->0) {
            crank.setGroup( pol[ rand() % pol.size() ] );
            sys+=crank.move(); // crank shaft move
          }
          break;
        case 5:
          k=pol.size();
          while (k-->0) {
            pivot.setGroup( pol[ rand() % pol.size() ] );
            sys+=pivot.move(); // pivot move
          }
          break;
      }

      for (auto i=pol.begin(); i!=pol.end()-1; i++)
        for (auto j=i+1; j!=pol.end(); j++)
          rdf( spc.geo->dist(i->cm,j->cm) )++;              // polymer mass-center distribution function

    } // end of micro loop

    sys.checkDrift( Energy::systemEnergy(spc,pot,spc.p)  ); // compare energy sum with current total

    cout << loop.timing();

  } // end of macro loop

  // save to disk
  rdf.save("rdf_p2p.dat");
  pqr.save("confout.pqr", spc.p);
  spc.save("state");

  // perform unit tests
  iso.test(test);
  gmv.test(test);
  mv.test(test);
  sys.test(test);

  // print information
  cout << loop.info() << sys.info() << mv.info() << gmv.info() << crank.info()
    << pivot.info() << iso.info() << shape.info() << test.info();
}
