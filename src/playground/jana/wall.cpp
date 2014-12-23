#include <faunus/faunus.h>

using namespace Faunus;
using namespace Faunus::Potential;

typedef Space<Geometry::Cuboidslit> Tspace;
//typedef CombinedPairPotential<DebyeHuckel,LennardJones> Tpairpot;
typedef CombinedPairPotential<DebyeHuckelShift,CutShift<LennardJones> > Tpairpot;

int main(int argc, char** argv) {

  InputMap mcp("wall.input");
  FormatXTC xtc(1000);                 // XTC gromacs trajectory format
  EnergyDrift sys;                     // class for tracking system energy drifts
  bool inPlane = mcp.get<bool>("molecule_plane");

  Tspace spc(mcp);

  // Surface-particle interactions
  auto surfpot = Energy::ExternalPotential<Tspace,Potential::GouyChapman<> >(mcp)
    + Energy::ExternalPotential<Tspace,Potential::HydrophobicWallLinear<> >(mcp);

  surfpot.first.expot.setSurfPositionZ( &spc.geo.len_half.z() );
  surfpot.second.expot.setSurfPositionZ( &spc.geo.len_half.z() );

  // Final energy function
  auto pot = Energy::NonbondedCutg2g<Tspace,Tpairpot>(mcp)
    + surfpot
    + Energy::EquilibriumEnergy<Tspace>(mcp);// + myenergy<Tspace>(mcp);

  //auto gouy = &pot.first.second;
  //gouy->expot.setSurfPositionZ( &spc.geo.len_half.z() ); // Pass position of GC surface

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
    FormatAAM::load(file,v);    // load molecule into 'v'
    Geometry::FindSpace fs;     // class for finding free space
    fs.allowMatterOverlap=true; // molecules can initially overlap
    if (inPlane)
      fs.dir=Point(1,1,0);      // molecule is in Z=0 plane
    fs.find(spc.geo,spc.p,v);   // place molecule randomly - coords saved in 'v'
    auto offset = Point(0,0,spc.geo.len_half.z()-mcp.get("molecule_offset",0));
    Geometry::translate(spc.geo, v, offset);
    pol[i] = spc.insert(v);
    pol[i].name=file;
    spc.enroll( pol[i] );
  }
  Group allpol( pol.front().front(), pol.back().back() );

  // Use charge scaling?
  bool boolChargeScaling=mcp.get<bool>("chargescaling",false);
  if (boolChargeScaling==true) {
    double D=pot.first.first.pairpot.first.debyeLength();
    for (auto &g : pol)
      for (auto i : g) {
        double kR = spc.p[i].radius/D;
        spc.p[i].charge *= std::sinh(kR)/kR;
        spc.trial[i] = spc.p[i];
      }
  }

  // Add atomic species
  Group salt;
  salt.addParticles(spc, mcp);
  salt.name="Free ions";
  spc.enroll(salt);

  Analysis::ChargeMultipole mpol;
  Average<double> P_ex;

  spc.load("state"); // load previous state, if any

  Move::TranslateRotate<Tspace> gmv(mcp,pot,spc);
  Move::AtomicTranslation<Tspace> mv(mcp,pot,spc);
  Move::SwapMove<Tspace> tit(mcp,pot,spc,*eqenergy);

  // MC translation only in XY plane?
  if (inPlane)
    for (auto &m : pol)
      gmv.directions[ m.name ]=Point(1,1,0);

  sys.init( Energy::systemEnergy(spc,pot,spc.p) );    // Store total system energy

  std::ofstream energyfile; // file for system energy
  energyfile.open("energy.dat");

  cout << atom.info() + spc.info() + pot.info() + tit.info()
    + textio::header("MC Simulation Begins!");

  int _numHydrophobic=0;
  for (auto &i : spc.p)
    if (i.hydrophobic)
      _numHydrophobic++;

  cout << "  Hydrophobic particles = " << _numHydrophobic << "\n"
    << "  Area per protein      = " << spc.geo.len.x()*spc.geo.len.y()/pol.size()
    << textio::_angstrom + textio::squared << "\n"
    << "  sqrt(area)            = " << sqrt(spc.geo.len.x()*spc.geo.len.y()/pol.size())
    << textio::_angstrom << "\n\n";

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

          if (slump()>0.999) {
            if (energyfile)
              energyfile << loop.count() << " " << sys.current()
                << " " << std::cbrt(spc.geo.getVolume()) << "\n";
          }
          break;
        case 2: // titration move
          sys+=tit.move();
          mpol.sample(pol,spc);
          break;
        case 3: // move free ions (if any)
          mv.setGroup(salt);
          sys+=mv.move();
          break;
      }

      if (slump()>0.99 ) {
        // 2D pressure calculation (see i.e. Frenkel+Smit)
        double sum=0;
        for (size_t i=0; i<pol.size()-1; i++) // double loop over proteins
          for (size_t j=i+1; j<pol.size(); j++)
            for (auto k : pol[i]) // loop over particles in proteins
              for (auto l : pol[j]) {
                auto r=spc.geo.vdist(spc.p[k], spc.p[l]); // distance vector
                auto f=pot.f_p2p(spc.p[k], spc.p[l]);     // force vector
                r.z()=f.z()=0;                            // neglect z-direction
                sum+=f.dot(r);                            // sum scalar product
              } 
        P_ex += sum / (2*spc.geo.len.x() * spc.geo.len.y());
      }

      if (slump()>0.999 ) {
        xtc.setbox(spc.geo.len);
        xtc.save("traj.xtc", spc.p); // save gromacs trajectory
      }
    } // end of micro loop

    sys.checkDrift( Energy::systemEnergy(spc,pot,spc.p) ); // detect energy drift

    cout << loop.timing();

  } // end of macro loop

  cout << tit.info() + loop.info() + sys.info() + gmv.info() + mv.info()
    + mpol.info();

  {
    using namespace textio;
    string punit=" kT/" + angstrom + squared;
    double P_id = pol.size() / (spc.geo.len.x()*spc.geo.len.y());
    cout << header("Analysis: 2D Pressure")
      << "  Excess              = " << P_ex.avg() << punit + "\n"
      << "  Ideal               = " << P_id << punit+"\n"
      << "  Total               = " << P_id+P_ex.avg() << punit+"\n"
      << "  Osmotic coefficient = " << 1 + P_ex.avg() / P_id << "\n\n";
  }

  FormatPQR::save("confout.pqr", spc.p, spc.geo.len);
  spc.save("state");
  mcp.save("mdout.mdp");
}
