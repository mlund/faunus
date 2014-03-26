#include <faunus/faunus.h>

/*
 * This will simulate:
 * - two rigid molecules fixed on the axis of
 *   a cylinder with open end planes and hard surfaces
 * - any number of atomic species
 * - equilibrium swap moves
 * - external pressure (NPT)
 * - electric multipole expansion analysis
 *
 * Particles interact with a combined Debye-Huckel/Lennard-Jones
 * potential with Lorentz-Berthelot mixed rules. Hydrophobic
 * groups interact with a custom epsilon as specified in the
 * input file.
 *
 * The structure of the molecules may be given either
 * as a PQR or AAM file -- see `FormatPQR` and `FormatAAM`
 * for details
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

  bool movie = mcp.get("movie", true);

  Tspace spc(mcp);
  auto pot = Energy::Nonbonded<Tspace,Tpairpot>(mcp)
    + Energy::ExternalPressure<Tspace>(mcp)
    + Energy::HydrophobicSASA<Tspace>(mcp)
    + Energy::EquilibriumEnergy<Tspace>(mcp);

  auto eqenergy = &pot.second;

  // set hydrophobic-hydrophobic LJ epsilon
  double epsh = mcp.get("eps_hydrophobic", 0.05);
  for (size_t i=0; i<atom.list.size()-1; i++)
    for (size_t j=i+1; j<atom.list.size(); j++)
      if (atom[i].hydrophobic)
        if (atom[j].hydrophobic)
          pot.first.first.first.pairpot.second.customEpsilon(i,j,epsh);

  // Add molecules
  int N1 = mcp.get("molecule1_N",0);
  int N2 = mcp.get("molecule2_N",0);
  bool inPlane = mcp.get("molecule_plane", false);
  string file;
  vector<Group> pol(N1+N2);
  for (int i=0; i<N1+N2; i++) {
    if (i>=N1)
      file = mcp.get<string>("molecule2");
    else
      file = mcp.get<string>("molecule1");
    Tspace::ParticleVector v;

    // PQR or AAM molecular file format?
    if (file.find(".pqr")!=std::string::npos)
      FormatPQR::load(file,v);
    else
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
  Analysis::MultipoleDistribution<Tspace> mpoldist(mcp);

  spc.load("state");

  Move::Isobaric<Tspace> iso(mcp,pot,spc);
  Move::TranslateRotateGroupCluster<Tspace> gmv(mcp,pot,spc);
  Move::AtomicTranslation<Tspace> mv(mcp, pot, spc);
  Move::SwapMove<Tspace> tit(mcp,pot,spc,*eqenergy);
  if (inPlane)
    for (auto &m : pol)
      gmv.directions[ m.name ]=Point(0,0,1);

  gmv.setMobile(salt);

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
          break;
        case 2:
          sys+=iso.move();
          break;
        case 3:
          mv.setGroup(salt);
          sys+=mv.move();
          break;
      }

      double xi = slp_global(); // random number [0,1)

      // Sample multipolar moments and distribution
      if ( xi>0.95 ) {
        pol[0].setMassCenter(spc);
        pol[1].setMassCenter(spc);
        mpol.sample(pol,spc);
        mpoldist.sample(spc, pol[0], pol[1]);
      }

      if ( movie==true && xi>0.995 ) {
        xtc.setbox( 1000. );
        xtc.save("traj.xtc", spc.p);
      }
    } // end of micro loop

    sys.checkDrift( Energy::systemEnergy(spc,pot,spc.p) ); // detect energy drift
    cout << loop.timing();

    spc.save("state");
    mcp.save("mdout.mdp");
    rdf.save("rdf_p2p.dat");
    mpoldist.save("multipole.dat");
    FormatPQR::save("confout.pqr", spc.p);

  } // end of macro loop

  cout << loop.info() + sys.info() + gmv.info() + mv.info()
    + iso.info() + tit.info() + mpol.info();

}
