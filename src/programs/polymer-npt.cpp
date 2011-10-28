#include <faunus/faunus.h>

using namespace Faunus;

typedef Geometry::Cuboid Tgeometry;    // select simulation geometry
typedef Potential::CoulombSR<Tgeometry, Potential::Coulomb, Potential::HardSphere> Tpairpot;

int main() {
  cout << textio::splash();
  atom.includefile("atomlist.inp");    // load atom properties
  InputMap mcp("polymer-npt.input");
  MCLoop loop(mcp);                    // class for handling mc loops
  FormatPQR pqr;                       // PQR structure file I/O
  FormatAAM aam;                       // AAM structure file I/O
  FormatXTC xtc(1000);                 // XTC gromacs trajectory format
  EnergyDrift sys;                     // class for tracking system energy drifts

  Energy::Hamiltonian pot;
  auto nonbonded = pot.create( Energy::Nonbonded<Tpairpot>(mcp) );
  auto bonded    = pot.create( Energy::Bonded() );
  Space spc( pot.getGeometry() );

  Move::Isobaric iso(mcp,pot,spc);
  Move::RotateGroup gmv(mcp,pot,spc);
  Move::ParticleTranslation mv(mcp, pot, spc);
  Analysis::PolymerShape shape;

  // Add salt
  GroupAtomic salt(spc, mcp);
  salt.name="Salt";

  // Add polymers
  vector<GroupMolecular> pol( mcp.get("polymer_N",0));
  string polyfile = mcp.get<string>("polymer_file", "");
  double req = mcp.get<double>("polymer_eqdist", 0);
  double k   = mcp.get<double>("polymer_forceconst", 0);
  Potential::Harmonic harmonic(k, req);
  atom["MM"].dp=10.;
  int i=1;
  for (auto &g : pol) {                    // load polymers
    aam.load(polyfile);
    g = spc.insert( aam.p );               // insert into space
    std::ostringstream o;
    o << "Polymer" << i++;
    g.name=o.str();
    spc.enroll(g);
    for (int i=g.beg; i<g.end; i++)
      bonded->bonds.add(i, i+1, harmonic); // add bonds
  }
  Group allpol( pol.front().beg, pol.back().end   );

  spc.load("space.state");

  double utot=pot.external();
  utot += pot.g_internal(spc.p, salt) + pot.g_external(spc.p, salt)
    + pot.g2g(spc.p, salt, allpol) + pot.g_internal(spc.p, allpol);
  for (auto &g : pol)
    utot += pot.g_external(spc.p, g);
  sys.init( utot );

  cout << atom.info() << spc.info() << pot.info() << mv.info()
    << textio::header("MC Simulation Begins!");

  while ( loop.macroCnt() ) {  // Markov chain 
    while ( loop.microCnt() ) {
      int k,i=rand() % 4;
      switch (i) {
        case 0:
          mv.setGroup(salt);
          sys+=mv.move();
          break;
        case 1:
          mv.setGroup(allpol);
          sys+=mv.move();
          for (auto &g : pol) {
            g.setMassCenter(spc);
            shape.sample(g,spc);
          }
          break;
        case 2:
          k=pol.size();
          while (k-->0) {
            gmv.setGroup( pol[ rand() % pol.size() ] );
            sys+=gmv.move();
          }
          break;
        case 3:
          sys+=iso.move();
          break;
      }
      if ( slp_global.runtest(0.01) ) // 1% chance
        xtc.save("traj.xtc", spc);

    } // end of micro loop

    double utot=pot.external();
    utot += pot.g_internal(spc.p, salt) + pot.g_external(spc.p, salt)
      + pot.g2g(spc.p, salt, allpol) + pot.g_internal(spc.p, allpol);
    for (auto &g : pol)
      utot += pot.g_external(spc.p, g);
    sys.checkDrift( utot );

    cout << loop.timing();
  } // end of macro loop

  pqr.save("confout.pqr", spc.p);
  spc.save("space.state");

  cout << loop.info() << sys.info() << mv.info() << gmv.info() << iso.info() << shape.info();
}
