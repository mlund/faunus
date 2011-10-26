#include <faunus/faunus.h>

using namespace Faunus;

typedef Geometry::Cuboid Tgeometry;    // select simulation geometry
typedef Potential::CoulombSR<Tgeometry, Potential::Coulomb, Potential::LennardJones> Tpairpot;

int main() {
  cout << textio::splash();
  atom.includefile("atomlist.inp");    // load atom properties
  InputMap mcp("bulk.inp");
  MCLoop loop(mcp);                    // class for handling mc loops
  FormatPQR pqr;                       // PQR structure file I/O
  FormatAAM aam;                       // AAM structure file I/O
  EnergyDrift sys;                     // class for tracking system energy drifts

  Energy::Hamiltonian pot;
  auto nonbonded = pot.create( Energy::Nonbonded<Tpairpot>(mcp) );
  auto bonded    = pot.create( Energy::Bonded() );
  Space spc( pot.getGeometry() );

  Move::Isobaric iso(mcp,pot,spc);
  Move::RotateGroup gmv(mcp,pot,spc);
  Move::ParticleTranslation mv(mcp, pot, spc);

  // Add particles
  GroupAtomic salt(spc, mcp);
  salt.name="Salt";
  mv.setGroup(salt);

  // Add polymers
  vector<Group> pol( mcp.get("Npolymers",0));
  for (auto &g : pol) {
    g = spc.insert( aam.p );
  }

  GroupMolecular test;
  test.beg=spc.p.size();
  spc.insert("Na", 2);
  test.end=spc.p.size()-1;
  test.name="Test";
  spc.p[test.beg].x=0;
  spc.p[test.beg].y=0;
  spc.p[test.beg].z=0;
  spc.p[test.end].x=2;
  spc.p[test.end].y=2;
  spc.p[test.end].z=2;
  spc.trial=spc.p;
  spc.enroll(test);

  spc.load("space.state");
  salt.setMassCenter(spc);
  test.setMassCenter(spc);

  //bonded->bonds.add(0,1, Potential::Harmonic(0.2, 10.0));
  //bonded->bonds.add(1,2, Potential::Harmonic(0.3,  5.0));

  #define UTOTAL \
      pot.g_external(spc.p, test)\
    + pot.g_internal(spc.p, salt)  + pot.g_external(spc.p, salt)\
    + pot.g2g(spc.p,salt,test) + pot.external()

  sys.init( UTOTAL );

  cout << atom.info() << spc.info() << pot.info() << mv.info()
       << textio::header("MC Simulation Begins!");

  while ( loop.macroCnt() ) {  // Markov chain 
    while ( loop.microCnt() ) {
      int i=rand() % 3;
      switch (i) {
        case 0:
          sys+=mv.move();
          break;
        case 1:
          gmv.setGroup(test);
          sys+=gmv.move();
          break;
        case 2:
          sys+=iso.move();
          break;
      }
    } // end of micro loop
    sys.checkDrift( UTOTAL );
    cout << loop.timing();
  } // end of macro loop

  pqr.save("confout.pqr", spc.p);
  spc.save("space.state");

  cout << loop.info() << sys.info() << mv.info() << gmv.info() << iso.info();
}
