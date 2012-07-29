#include <faunus/faunus.h>

using namespace Faunus;

typedef Geometry::Sphere Tgeometry;
typedef Potential::DebyeHuckelLJ Tpairpot;

int main() {
  InputMap mcp("input");
  MCLoop loop(mcp);                    // class for handling mc loops
  FormatPQR pqr;                       // PQR structure file I/O
  FormatAAM aam;                       // AAM structure file I/O
  EnergyDrift sys;                     // class for tracking system energy drifts
  UnitTest test(mcp);

  Energy::Hamiltonian pot;
  auto nonbonded = pot.create( Energy::Nonbonded<Tpairpot,Tgeometry>(mcp) );
  Space spc( pot.getGeometry() );

  Move::SwapMove tit(mcp,pot,spc);
  Analysis::ChargeMultipole poleTotal;
  
  // Add molecule to middle of simulation container
  string file = mcp.get<string>("molecule", "in.aam");
  aam.load(file);
  Point cm = Geometry::massCenter(*spc.geo, aam.particles() );
  Geometry::translate(*spc.geo, aam.particles(), -cm); // place in origo
  GroupMolecular g = spc.insert(aam.particles());
  g.name=file;
  spc.enroll(g);

  for(unsigned int i=0; i < spc.p.size(); i++){
    spc.p[i].charge=atom[spc.p[i].id].charge;
    spc.trial[i].charge=spc.p[i].charge;
  }
  
  tit.findSites(spc.p);  // search for titratable sites
  spc.load("state");

  double utot=pot.external() + pot.g_internal(spc.p, g);
  sys.init( utot );

  cout << atom.info() << spc.info() << pot.info() << tit.info()
    << textio::header("MC Simulation Begins!");

  while ( loop.macroCnt() ) {  // Markov chain 
    while ( loop.microCnt() ) {
      int i=rand() % 1;
      switch (i) {
        case 0:
          sys+=tit.move();
          poleTotal.sample(g,spc);
          break;
      }
    } // end of micro loop

    double utot=pot.external() + pot.g_internal(spc.p, g);
    sys.checkDrift( utot );
    cout << loop.timing() << std::flush;

  } // end of macro loop

  cout << loop.info() << sys.info() << tit.info() << g.info() << endl 
       << textio::header("Total Charge Analysis") << poleTotal.info();

  pqr.save("confout.pqr", spc.p);
  spc.save("state");
}
