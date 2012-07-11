/*
 * manybodyMPI.cpp
 *
 * This will simulate a mix of rigid molecules and atomic species in a cubic simulation
 * container. The default Hamiltonian is:
 * - Debye-Huckel electrostatics
 * - Lennard-Jones interactions
 * - Gouy-Chapman electrostatics (for slit containers)
 *
 * The following Monte Carlo Moves are implemented:
 * 1) Translation of atomic species and rigid molecules
 * 2) Rotation of rigid molecules
 * 3) Particle swap moves (proton titration etc.)
 * 4) Isobaric volume moves (NPT ensemble)
 * 5) Parallel Tempering using MPI
 */

#include <faunus/faunus.h>

using namespace Faunus;

#ifdef SLIT
typedef Geometry::Cuboidslit Tgeometry;
#else
typedef Geometry::Cuboid Tgeometry;
#endif
typedef Potential::CombinedPairPotential<Potential::DebyeHuckel, Potential::LennardJones> Tpairpot;

int main(int argc, char** argv) {
  Faunus::MPI::MPIController mpi;
  mpi.cout << textio::splash();

  InputMap mcp(textio::prefix+"manybody.input");
  MCLoop loop(mcp);                    // class for handling mc loops
  FormatPQR pqr;                       // PQR structure file I/O
  FormatAAM aam;                       // AAM structure file I/O
  FormatTopology top;
  FormatXTC xtc(1000);                 // XTC gromacs trajectory format
  EnergyDrift sys;                     // class for tracking system energy drifts
  UnitTest test(mcp);

  Energy::Hamiltonian pot;
  auto nonbonded = pot.create( Energy::Nonbonded<Tpairpot,Tgeometry>(mcp) );
#ifdef SLIT
  auto gouy = pot.create( Energy::GouyChapman(mcp) );
  gouy->setPosition( nonbonded->geometry.len_half.z );
#endif
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
    pol[i].name=file;
    spc.enroll( pol[i] );
  }
  Group allpol( pol.front().front(), pol.back().back() );

  // Add atomic species
  GroupAtomic salt(spc, mcp);
  salt.name="Atomic Species";
  spc.enroll(salt);

  Move::Isobaric iso(mcp,pot,spc);
  Move::TranslateRotate gmv(mcp,pot,spc);
  Move::AtomicTranslation mv(mcp, pot, spc);
  Move::SwapMove tit(mcp,pot,spc);
  Move::ParallelTempering temper(mcp,pot,spc,mpi);

  Analysis::RadialDistribution<float,int> rdf(0.25);
  Analysis::ChargeMultipole mpol;

  spc.load(textio::prefix+"state");

  sys.init( Energy::systemEnergy(spc,pot,spc.p) );

  mpi.cout << atom.info() << spc.info() << pot.info() << tit.info()
    << textio::header("MC Simulation Begins!");

  while ( loop.macroCnt() ) {  // Markov chain 
    while ( loop.microCnt() ) {
      int k,i=rand() % 4;
      switch (i) {
        case 0: // translate and rotate molecules
          k=pol.size();
          while (k-->0) {
            gmv.setGroup( pol[ rand() % pol.size() ] );
            sys+=gmv.move();
          }
          for (auto i=pol.begin(); i!=pol.end()-1; i++)
            for (auto j=i+1; j!=pol.end(); j++)
              rdf( spc.geo->dist(i->cm,j->cm) )++;
          break;
        case 1: // volume move
          sys+=iso.move();
          break;
        case 2: // titration move
          sys+=tit.move();
          mpol.sample(pol,spc);
          break;
        case 3: // translate atomic species
          mv.setGroup(salt);
          sys+=mv.move();
          break;
      }

      if ( slp_global.runtest(0.0001) ) {
        xtc.setbox( nonbonded->geometry.len );
        xtc.save(textio::prefix+"traj.xtc", spc);
      }
    } // end of micro loop

    temper.setCurrentEnergy( sys.current() );
    sys+=temper.move();

    sys.checkDrift( Energy::systemEnergy(spc,pot,spc.p) );

    mpi.cout << loop.timing();

  } // end of macro loop

  //iso.test(test);
  //gmv.test(test);
  //sys.test(test);

  mpi.cout << loop.info() << sys.info() << gmv.info() << mv.info() << iso.info() << tit.info()
    << mpol.info() << temper.info();

  rdf.save(textio::prefix+"rdf_p2p.dat");
  pqr.save(textio::prefix+"confout.pqr", spc.p);
  //top.save(textio::prefix+"mytopol.top", spc);
  mcp.save(textio::prefix+"mdout.mdp");
  spc.save(textio::prefix+"state");
}
