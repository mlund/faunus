#include <faunus/faunus.h>
//#include "range.hpp"

using namespace Faunus;

typedef Geometry::Cuboid Tgeometry;                // select simulation geometry
typedef Potential::CoulombHS Tpairpot;

int main() {
  Group g(100,1);
  g.resize(0);
  cout << "size = " << g.size() << endl;
  cout << "range = " << g.front() << " " << g.back() << endl;
  cout << "empty bool = " << g.empty() << endl;
  for (auto i=g.begin(); i!=g.end(); i++)
    cout << *i << endl;

  cout << textio::splash();
  atom.includefile("atomlist.inp");    // load atom properties
  InputMap mcp("bulk.inp");
  MCLoop loop(mcp);                    // class for handling mc loops
  FormatPQR pqr;                       // PQR structure file I/O
  EnergyDrift sys;                     // class for tracking system energy drifts

  Energy::Hamiltonian pot;
  auto nonbonded = pot.create( Energy::Nonbonded<Tpairpot,Tgeometry>(mcp) );
  Space spc( pot.getGeometry() );

  // Handle particles
  GroupAtomic salt(spc, mcp);
  salt.name="Salt";
  spc.load("space.state", Space::RESIZE);

  Move::GrandCanonicalSalt gc(mcp,pot,spc,salt);
  Move::AtomicTranslation mv(mcp, pot, spc);  // Particle move class
  mv.setGroup(salt);

  // Widom particle insertion
  Analysis::RadialDistribution<float,unsigned int> rdf(0.2);
  rdf.maxdist = 50.;
  particle a;
  Analysis::Widom widom(spc, pot);
  a = atom["Cl"];
  widom.addGhost( spc);
  widom.addGhost( a );

#define UTOTAL \
  + pot.g_internal(spc.p, salt)  + pot.g_external(spc.p, salt)\
  + pot.external()
  sys.init( UTOTAL );

  cout << atom.info() << spc.info() << pot.info() << mv.info()
    << textio::header("MC Simulation Begins!");

  while ( loop.macroCnt() ) {  // Markov chain 
    while ( loop.microCnt() ) {
      sys+=mv.move( salt.size() );
      sys+=gc.move( salt.size()/2 );
      if (slp_global.randOne()>0.9) {
        widom.sample();
        rdf.sample(spc,salt,atom["Na"].id, atom["Cl"].id);
      }
    }
    sys.checkDrift( UTOTAL );
    cout << loop.timing();
  }

  pqr.save("confout.pqr", spc.p);
  spc.save("space.state");
  rdf.save("rdf.dat");

  cout << sys.info() << loop.info() << mv.info() << gc.info() << widom.info();
}
