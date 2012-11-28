#include <faunus/faunus.h>
#include <tclap/CmdLine.h>

using namespace Faunus;

#define CYLINDER

#ifdef CYLINDER
typedef Geometry::PeriodicCylinder Tgeometry;
#else
typedef Geometry::Cuboid Tgeometry;
#endif
typedef Potential::CoulombSR<Tgeometry, Potential::DebyeHuckel, Potential::SoftRepulsion> Tpairpot;

int main(int argc, char** argv) {
  cout << textio::splash();
  InputMap mcp("cluster.input");
  MCLoop loop(mcp);                    // class for handling mc loops
  FormatPQR pqr;                       // PQR structure file I/O
  FormatAAM aam;                       // AAM structure file I/O
  FormatTopology top;
  FormatXTC xtc(1000);                 // XTC gromacs trajectory format
  EnergyDrift sys;                     // class for tracking system energy drifts
  UnitTest test(mcp);

  Energy::Hamiltonian pot;
  auto nonbonded = pot.create( Energy::Nonbonded<Tpairpot>(mcp) );
  auto constrain = pot.create( Energy::MassCenterConstrain(pot.getGeometry()) );
  Space spc( pot.getGeometry() );

  // Add molecular species
  vector<GroupMolecular> pol( mcp.get("polymer_N",0));
  string polyfile = mcp.get<string>("polymer_file", "");
  for (auto &g : pol) {
    aam.load(polyfile);
    Geometry::FindSpace f;
    f.dir.x()=0; // put mass center
    f.dir.y()=0; //   at [x,y,z] = [0,0,random]
    if (f.find(*spc.geo, spc.p, aam.p )) {
      g = spc.insert( aam.p );
      g.name="Protein";
      spc.enroll(g);
    } else
      return 1;
  }

  constrain->addPair(pol[0], pol[1], 20, 80);

  // Add salt
  GroupAtomic salt(spc, mcp);
  salt.name="Salt";
  spc.enroll(salt);
  spc.load("state");

  Move::TranslateRotateCluster gmv(mcp,pot,spc);
  Move::AtomicTranslation mv(mcp, pot, spc);
  gmv.setMobile(salt); // specify where to look for clustered ions
  mv.setGroup(salt);   // specify atomic particles to be moved

#ifdef CYLINDER
  gmv.dir.x()=0; // do not move in x
  gmv.dir.y()=0; // nor y direction
  Analysis::LineDistribution<float,unsigned long int> rdf(0.25);
#else
  Analysis::RadialDistribution<float,int> rdf(0.25);
#endif
  rdf.maxdist=80;

  double utot=pot.external() + pot.g_internal(spc.p, salt);
  for (auto &p : pol)
    utot+=pot.g2g(spc.p, salt, p);
  for (auto i=pol.begin(); i!=pol.end()-1; i++)
    for (auto j=i+1; j!=pol.end(); j++)
      utot += pot.g2g(spc.p, *i, *j);
  sys.init( utot );

  cout << atom.info() << spc.info() << pot.info() << textio::header("MC Simulation Begins!");

  while ( loop.macroCnt() ) {  // Markov chain 
    while ( loop.microCnt() ) {
      int k,i=rand() % 2;
      switch (i) {
        case 0:
          mv.setGroup(salt);
          sys+=mv.move( salt.size()/2 );
          break;
        case 1:
          k=pol.size();
          while (k-->0) {
            gmv.setGroup( pol[ rand() % pol.size() ] );
            sys+=gmv.move();
          }
          for (auto i=pol.begin(); i!=pol.end()-1; i++)
            for (auto j=i+1; j!=pol.end(); j++) {
              double r=spc.geo->dist(i->cm,j->cm);
              if (r<rdf.maxdist) {
                rdf(r)++;
              }
            }
          break;
      }
    } // end of micro loop

    double utot=pot.external() + pot.g_internal(spc.p, salt);
    for (auto &p : pol)
      utot+=pot.g2g(spc.p, salt, p);
    for (auto i=pol.begin(); i!=pol.end()-1; i++)
      for (auto j=i+1; j!=pol.end(); j++)
        utot += pot.g2g(spc.p, *i, *j);
    sys.checkDrift( utot );
    rdf.save("rdf_p2p.dat");
    cout << loop.timing();
  } // end of macro loop

  pqr.save("confout.pqr", spc.p);
  spc.save("state");

  cout << loop.info() << spc.info() << sys.info() << mv.info() << gmv.info();
}
