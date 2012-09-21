#include <faunus/faunus.h>

using namespace Faunus;
using namespace Faunus::Potential;

typedef Geometry::Cuboid Tgeometry;
typedef CombinedPairPotentialCutoff<DebyeHuckelShift,LennardJones> Tpairpot;
typedef CombinedPairPotentialCutoff<DebyeHuckel,LennardJones> Tpairpot2;
//typedef Potential::DebyeHuckelLJ Tpairpot;

int main(int argc, char** argv) {
  string inputfile="manybody.input";
  string state="state";

  InputMap mcp(inputfile);
  MCLoop loop(mcp);                    // class for handling mc loops
  FormatPQR pqr;                       // PQR structure file I/O
  FormatGRO gro;
  FormatAAM aam;                       // AAM structure file I/O
  FormatTopology top;
  FormatXTC xtc(1000);                 // XTC gromacs trajectory format
  EnergyDrift sys;                     // class for tracking system energy drifts
  UnitTest test(mcp);

  Energy::Hamiltonian pot;
  auto nonbonded = pot.create( Energy::Nonbonded<Tpairpot,Tgeometry>(mcp) );
  Space spc( pot.getGeometry() );

  // Add molecules
  int N1 = mcp.get("molecule_N1",0);
  int N2 = mcp.get("molecule_N2",0);
  string file;
  vector<GroupMolecular> pol(N1+N2);
  for (int i=0; i<N1+N2; i++) {
    GroupMolecular g;
    if (i>=N1)
      file = mcp.get<string>("molecule_file2", "");
    else
      file = mcp.get<string>("molecule_file1", "");
    aam.load(file);
    Geometry::FindSpace f;
    f.find(*spc.geo, spc.p, aam.p);        // find empty spot in particle vector
    pol[i] = spc.insert( aam.p );          // insert into space
    pol[i].name="prot";
    spc.enroll( pol[i] );
  }
  Group allpol( pol.front().front(), pol.back().back() );

  spc.load(state);

  Move::Isobaric iso(mcp,pot,spc);
  Move::TranslateRotate gmv(mcp,pot,spc);
  Move::SwapMove tit(mcp,pot,spc);

  Analysis::RadialDistribution<float,int> rdf(0.25);
  Analysis::ChargeMultipole mpol;

  sys.init( Energy::systemEnergy(spc,pot,spc.p) );

  cout << atom.info() << spc.info() << pot.info() << tit.info()
    << textio::header("MC Simulation Begins!");

  while ( loop.macroCnt() ) {  // Markov chain 
    while ( loop.microCnt() ) {
      int k,i=rand() % 3;
      switch (i) {
        case 0:
          k=pol.size();
          while (k-->0) {
            gmv.setGroup( pol[ rand() % pol.size() ] );
            sys+=gmv.move();
          }
          for (auto i=pol.begin(); i!=pol.end()-1; i++)
            for (auto j=i+1; j!=pol.end(); j++)
              rdf( spc.geo->dist(i->cm,j->cm) )++;
          break;
        case 1:
          sys+=iso.move();
          break;
        case 2:
          sys+=tit.move();
          mpol.sample(pol,spc);
          break;
      }
      if ( slp_global.runtest(0.001) ) {
        xtc.setbox( nonbonded->geometry.len );
        xtc.save("traj.xtc", spc.p);
      }
    } // end of micro loop

    sys.checkDrift( Energy::systemEnergy(spc,pot,spc.p) );

    cout << loop.timing();

  } // end of macro loop

  /*
  iso.test(test);
  gmv.test(test);
  sys.test(test);
  */

  cout << loop.info() << sys.info() << gmv.info() << iso.info() << tit.info()
    << mpol.info();

  rdf.save("rdf_p2p.dat");
  pqr.save("confout.pqr", spc.p);
  gro.len = nonbonded->geometry.len.x;
  gro.save("confout.gro", spc);
  top.save("mytopol.top", spc);

  nonbonded->pairpot.save("pairpot-NaCl.dat", atom["Na"].id, atom["Cl"].id);
  Tpairpot2 pot2(mcp);
  pot2.save("pairpot-full-NaCl.dat", atom["Na"].id, atom["Cl"].id);

  spc.save(state);
}
