/*
 * manybodyMPI.cpp
 *
 * This will simulate polymers in a rectangular slit container:
 * container. The default Hamiltonian is:
 * - Debye-Huckel electrostatics
 * - Lennard-Jones interactions
 * - Gouy-Chapman electrostatics (for slit containers)
 *
 * The following Monte Carlo Moves are implemented:
 * 1) Translation / rotation of rigid molecules
 * 2) Particle swap moves (proton titration etc.)
 * 3) Isobaric volume moves (NPT ensemble)
 * 4) Parallel Tempering using MPI
 */

#include <faunus/faunus.h>

using namespace Faunus;

#ifdef SLIT
typedef Geometry::Cuboidslit Tgeometry;
#else
//typedef Geometry::Sphere Tgeometry;
typedef Geometry::Cuboid Tgeometry;
#endif
typedef Potential::CombinedPairPotential<Potential::Harmonic, Potential::LennardJones> Tbondpot;
typedef Potential::DebyeHuckelLJ Tpairpot;
//typedef Potential::DebyeHuckelr12 Tpairpot;

int main(int argc, char** argv) {
  Faunus::MPI::MPIController mpi;
  mpi.cout << textio::splash();

  InputMap mcp(textio::prefix+"manybody.input");
  MCLoop loop(mcp);                    // class for handling mc loops
  FormatPQR pqr;                       // PQR structure file I/O
  FormatAAM aam;                       // AAM structure file I/O
  FormatTopology top;
  FormatXTC xtc(1000);                 // XTC gromacs trajectory format for only cuboid geomtries
  EnergyDrift sys;                     // class for tracking system energy drifts
  UnitTest test(mcp);

  Energy::Hamiltonian pot;
  auto nonbonded = pot.create( Energy::Nonbonded<Tpairpot,Tgeometry>(mcp) );
  auto bonded    = pot.create( Energy::Bonded() );
#ifdef SLIT
  auto gouy = pot.create( Energy::GouyChapman(mcp) );
  gouy->setPosition( nonbonded->geometry.len_half.z );
#endif
  Space spc( pot.getGeometry() );

  // Add polymers
  vector<GroupMolecular> pol( mcp.get("polymer_N",0) );
  string polyfile = mcp.get<string>("polymer_file", "");
  atom["MM"].dp = 0.;
  int ii=1;
  mpi.cout << "Number of polymers: " << pol.size() << endl;
  for (auto &g : pol) {                    // load polymers
    aam.load(polyfile);
    Geometry::FindSpace f;
    f.find(*spc.geo, spc.p, aam.particles());        // find empty spot in particle vector
    g = spc.insert( aam.particles() );   // insert into space
    ostringstream o;
    o << "Polymer" << ii++;
    g.name=o.str();
    spc.enroll(g);
    for (int i=g.front(); i<g.back(); i++)
      bonded->add(i, i+1, Tbondpot(mcp, "polymer_", "minuslj_")); // add bonds
  }
  Group allpol( pol.front().front(), pol.back().back() ); // make group w. all polymers
  
  // Add salt
  GroupAtomic ball(spc, mcp);
  ball.name="ball";
  mpi.cout << "Number of balls: " << ball.size() << endl;

  // atom["NTR"].dp = 10.;
  // atom["CTR"].dp = 10.;
  // atom["HIS"].dp = 10.;
  // atom["HNTR"].dp = 10.;
  // atom["HCTR"].dp = 10.;
  // atom["HHIS"].dp = 10.;

  //Move::Isobaric iso(mcp,pot,spc);
  Move::TranslateRotate gmv(mcp,pot,spc);
  Move::AtomicTranslation mv(mcp, pot, spc);
  Move::AtomicTranslation mv2(mcp, pot, spc, "ball");
  Move::SwapMove tit(mcp,pot,spc);
  Move::CrankShaft crank(mcp, pot, spc);
  Move::Pivot pivot(mcp, pot, spc);
  Move::ParallelTempering temper(mcp,pot,spc,mpi);

  Analysis::RadialDistribution<float,int> rdf(0.25);
  Analysis::LineDistribution<float,int> surfdist(0.25);
  Analysis::ChargeMultipole mpol;
  Analysis::PolymerShape shape;
  std::map< string , Analysis::LineDistribution<float,int> > surfmap;

  spc.load(textio::prefix+"state");

  sys.init( Energy::systemEnergy(spc,pot,spc.p) );

  mpi.cout << atom.info() << spc.info() << pot.info() //<< tit.info()
    << textio::header("MC Simulation Begins!");

  pqr.save(textio::prefix+"initial.pqr", spc.p);
  while ( loop.macroCnt() ) {  // Markov chain 
    while ( loop.microCnt() ) {
      int k,i=rand() % 5;
      switch (i) {
        case 0: // translate and rotate molecules
          k=pol.size();//+ball.size();
          while (k-->0) {
            int s=rand()%2;
            if (s==1) {
              gmv.setGroup( pol.at( rand() % pol.size() ) );
              sys+=gmv.move(); 
            }
            else {
              mv2.setGroup(ball); //atomic translate move
              sys+=mv2.move(1); 
            }
          }
#ifdef SLIT
          for (auto &g : pol) {
            surfmap[g.name]( gouy->dist2surf(g.cm) )++;
            surfdist( gouy->dist2surf(g.cm) )++; // polymer mass center to GC surface histogram
          }
#endif
          for (auto i=pol.begin(); i!=pol.end()-1; i++)
            for (auto j=i+1; j!=pol.end(); j++)
              rdf( spc.geo->dist(i->cm,j->cm) )++; // sample polymer-polymer rdf
          break;
        case 1:
          mv.setGroup(allpol);
          sys+=mv.move( allpol.size() ); // translate monomers
          for (auto &g : pol) {
            g.setMassCenter(spc);
            shape.sample(g,spc);
          }
          break;
        case 2: // titration move
          sys+=tit.move();
          mpol.sample(pol,spc);
          break;
        case 3: // crankshaft move
          k=pol.size();
          while (k-->0) {
            crank.setGroup( pol[ rand() % pol.size() ] );
            sys+=crank.move();
          }
          for (auto &g : pol) {
            g.setMassCenter(spc);
            shape.sample(g,spc);
          }
          break;
        case 4: //pivot move
          k=pol.size();
          while (k-->0) {
            pivot.setGroup( pol[ rand() % pol.size() ] );
            sys+=pivot.move();
          }
          for (auto &g : pol) {
            g.setMassCenter(spc);
            shape.sample(g,spc);
          }
          break;
        case 5: // volume move
          //sys+=iso.move();
          break;
        }

      if ( slp_global.runtest(0.0001) ) {
        xtc.setbox( nonbonded->geometry.len );
        xtc.save(textio::prefix+"traj.xtc", spc);
      }
#ifdef SLIT
      if ( slp_global.runtest(0.0001) ) {
        

      }
#endif
    } // end of micro loop

    temper.setCurrentEnergy( sys.current() );
    sys+=temper.move();

    sys.checkDrift( Energy::systemEnergy(spc,pot,spc.p) );

    mpi.cout << loop.timing();
    //cout << loop.timing();

  } // end of macro loop

  mpi.cout << loop.info() << sys.info() << gmv.info() << tit.info() << mv.info() << mv2.info() 
    << crank.info() << pivot.info() << temper.info() //<< iso.info()
    << mpol.info()  << shape.info();   

  rdf.save(textio::prefix+"rdf_p2p.dat");
  surfdist.save(textio::prefix+"surfdist.dat");
  pqr.save(textio::prefix+"confout.pqr", spc.p);
  mcp.save(textio::prefix+"mdout.mdp");
  spc.save(textio::prefix+"state");
  for (auto &g : pol) 
    surfmap[g.name].save(textio::prefix+g.name+"surfdist.dat");
}
