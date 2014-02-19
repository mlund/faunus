#include <faunus/faunus.h>

/*
 * This will simulate:
 * - an arbitrary number of rigid molecules (two different kinds)
 * - any number of atomic species
 * - equilibrium swap moves
 * - external pressure (NPT)
 */

using namespace Faunus;
using namespace Faunus::Potential;

typedef Space<Geometry::PeriodicCylinder> Tspace;
typedef CombinedPairPotential<DebyeHuckel,LennardJonesLB> Tpairpot;

int main(int argc, char** argv) {
  InputMap mcp("cyl.input");
  MCLoop loop(mcp);                    // class for handling mc loops
  FormatXTC xtc(1000);                 // XTC gromacs trajectory format
  EnergyDrift sys;                     // class for tracking system energy drifts
  UnitTest test(mcp);
  bool inPlane = mcp.get<bool>("molecule_plane");

  Tspace spc(mcp);
  auto pot = Energy::Nonbonded<Tspace,Tpairpot>(mcp)
    + Energy::ExternalPressure<Tspace>(mcp)
    + Energy::EquilibriumEnergy<Tspace>(mcp);

  auto eqenergy = &pot.second;

  // Add molecules
  int N1 = mcp.get("molecule1_N",0);
  int N2 = mcp.get("molecule2_N",0);
  string file;
  vector<Group> pol(N1+N2);
  for (int i=0; i<N1+N2; i++) {
    if (i>=N1)
      file = mcp.get<string>("molecule2");
    else
      file = mcp.get<string>("molecule1");
    Tspace::ParticleVector v;
    FormatAAM::load(file,v);
    Geometry::FindSpace fs;
    if (inPlane)
      fs.dir=Point(0,0,1);
    fs.find(spc.geo,spc.p,v);
    pol[i] = spc.insert(v);
    pol[i].name=file+std::to_string(i);
    spc.enroll( pol[i] );
  }
  Group allpol( pol.front().front(), pol.back().back() );

  // Add atomic species
  Group salt;
  salt.addParticles(spc, mcp);
  salt.name="Atomic Species";
  spc.enroll(salt);

  Analysis::LineDistribution<> rdf(0.5);
  Analysis::ChargeMultipole mpol;

  spc.load("state");

  Move::Isobaric<Tspace> iso(mcp,pot,spc);
  Move::TranslateRotate<Tspace> gmv(mcp,pot,spc);
  Move::AtomicTranslation<Tspace> mv(mcp, pot, spc);
  Move::SwapMove<Tspace> tit(mcp,pot,spc,*eqenergy);
  if (inPlane)
    for (auto &m : pol)
      gmv.directions[ m.name ]=Point(0,0,1);

  eqenergy->eq.findSites(spc.p);

  sys.init( Energy::systemEnergy(spc,pot,spc.p) );    // Store total system energy

  cout << atom.info() << spc.info() << pot.info() << tit.info()
    << textio::header("MC Simulation Begins!");

  while ( loop.macroCnt() ) {  // Markov chain 
    while ( loop.microCnt() ) {
      int k,i=rand() % 2;
      switch (i) {
        case 0:
          k=pol.size();
          while (k-->0) {
            gmv.setGroup( pol[ rand() % pol.size() ] );
            sys+=gmv.move();
          }
          for (auto i=pol.begin(); i!=pol.end()-1; i++)
            for (auto j=i+1; j!=pol.end(); j++)
              rdf( spc.geo.dist(i->cm,j->cm) )++;
          break;
        case 1:
          sys+=tit.move();
          if( slp_global()>0.9 )
            mpol.sample(pol,spc);
          break;
        case 2:
          sys+=iso.move();
          break;
        case 3:
          mv.setGroup(salt);
          sys+=mv.move();
          break;
      }
      if ( slp_global()<-0.001 ) {
        xtc.setbox( 1000. );
        xtc.save("traj.xtc", spc);
      }
    } // end of micro loop

    sys.checkDrift( Energy::systemEnergy(spc,pot,spc.p) ); // detect energy drift
    cout << loop.timing();

  } // end of macro loop

  cout << loop.info() + sys.info() + gmv.info() + mv.info()
    + iso.info() + tit.info() + mpol.info();

  rdf.save("rdf_p2p.dat");
  FormatPQR::save("confout.pqr", spc.p);
  spc.save("state");
  mcp.save("mdout.mdp");
}
