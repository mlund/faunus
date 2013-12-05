#include <faunus/faunus.h>

using namespace Faunus;
using namespace Faunus::Potential;

typedef Space<Geometry::Cuboidslit> Tspace;
typedef CombinedPairPotential<DebyeHuckel,LennardJones> Tpairpot;
//typedef CombinedPairPotential<DebyeHuckelShift,CutShift<LennardJones> > Tpairpot;

int main(int argc, char** argv) {

  InputMap mcp("wall.input");
  FormatXTC xtc(1000);                 // XTC gromacs trajectory format
  EnergyDrift sys;                     // class for tracking system energy drifts
  bool inPlane = mcp.get<bool>("molecule_plane");

  Tspace spc(mcp);
  auto pot = Energy::Nonbonded<Tspace,Tpairpot>(mcp)
    + Energy::ExternalPotential<Tspace,Potential::GouyChapman<> >(mcp)
    + Energy::EquilibriumEnergy<Tspace>(mcp);// + myenergy<Tspace>(mcp);

  auto gouy = &pot.first.second;
  gouy->expot.setSurfPositionZ( &spc.geo.len_half.z() ); // Pass position of GC surface

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
    fs.allowMatterOverlap=true;
    if (inPlane)
      fs.dir=Point(1,1,0);
    fs.find(spc.geo,spc.p,v);
    pol[i] = spc.insert(v);
    pol[i].name=file;
    spc.enroll( pol[i] );
  }
  Group allpol( pol.front().front(), pol.back().back() );

  // Add atomic species
  Group salt;
  salt.addParticles(spc, mcp);
  salt.name="Atomic Species";
  spc.enroll(salt);

  Analysis::LineDistribution<> cmhist(0.2);            // monomer-surface histogram
  Analysis::ChargeMultipole mpol;

  spc.load("state"); // load previous state, if any

  Move::TranslateRotate<Tspace> gmv(mcp,pot,spc);
  Move::AtomicTranslation<Tspace> mv(mcp,pot,spc);
  Move::SwapMove<Tspace> tit(mcp,pot,spc,*eqenergy);

  if (inPlane)
    for (auto &m : pol)
      gmv.directions[ m.name ]=Point(1,1,0);

  sys.init( Energy::systemEnergy(spc,pot,spc.p) );    // Store total system energy

  std::ofstream energyfile;

  energyfile.open("energy.dat");
  cout << atom.info() << spc.info() << pot.info() << tit.info()
    << textio::header("MC Simulation Begins!");

  MCLoop loop(mcp);
  while ( loop.macroCnt() ) {  // Markov chain
    while ( loop.microCnt() ) {
      int k,i=rand() % 3;
      switch (i) {
        case 1: // move all proteins
          k=pol.size();
          while (k-->0) {
            gmv.setGroup( pol[ rand() % pol.size() ] );
            sys+=gmv.move();
          }

          cmhist(  gouy->expot.surfDist(pol[0].cm)  )++;

          if (slp_global()>0.99) {
            if (energyfile)
              energyfile << loop.count() << " " << sys.current()
                << " " << std::cbrt(spc.geo.getVolume()) << "\n";
          }
          break;
        case 2: // titration move
          sys+=tit.move();
          mpol.sample(pol,spc);
          break;
        case 3: // move atomic species
          mv.setGroup(salt);
          sys+=mv.move();
          break;
      }
      if (slp_global()>0.95 ) {
        xtc.setbox( spc.geo.len );
        xtc.save("traj.xtc", spc);
      }
    } // end of micro loop

    sys.checkDrift( Energy::systemEnergy(spc,pot,spc.p) ); // detect energy drift

    cout << loop.timing();

  } // end of macro loop

  cout << tit.info() + loop.info() + sys.info() + gmv.info() + mv.info()
    + mpol.info() << endl;

  cmhist.save("cmdist.dat");
  FormatPQR::save("confout.pqr", spc.p, spc.geo.len);
  spc.save("state");
  mcp.save("mdout.mdp");
}
