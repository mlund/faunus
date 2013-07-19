#ifdef ENABLE_MPI

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

template<class Tspace>
struct myenergy : public Energy::Energybase<Tspace> {

  typedef Energy::Energybase<Tspace> base;

  double rc,a,z;
  bool havecharge;

  Potential::DebyeHuckel dh;

  myenergy(InputMap &in) : dh(in) {
    base::name="myenergy";
    havecharge=false;
    double k=1/dh.debyeLength();
    rc=in.get<double>("dh_cutoff", pc::infty);
    double v=4*pc::pi*pow(rc,3)/3.; // cut-off volume
    a=4*pc::pi*exp(-k*rc)*(k*rc+1)/(k*k);
    a=a*dh.bjerrumLength() / v;
  }

  double i_external(const p_vec &p, int i) FOVERRIDE {
    if (havecharge==false) {
      z=base::spc->charge();
      havecharge=true;
    }
    return p[i].charge * z * a;
  }

  double g_external(const p_vec &p, Group &g) FOVERRIDE {
    double u=0;
    for (auto i : g)
      u+=i_external(p,i);
    return u;
  }

  string _info() {
    return "myenergy - cutoff = " + std::to_string(rc)
      + textio::_angstrom;
  }
};

typedef Space<Geometry::Cuboid> Tspace;
typedef CombinedPairPotential<DebyeHuckelShift,CutShift<LennardJones> > Tpairpot;

int main(int argc, char** argv) {
  Faunus::MPI::MPIController mpi;
  
  InputMap mcp("manybody.input");
  FormatXTC xtc(1000);                 // XTC gromacs trajectory format
  EnergyDrift sys;                     // class for tracking system energy drifts
  UnitTest test(mcp);
  bool inPlane = mcp.get<bool>("molecule_plane");

  Tspace spc(mcp);
  auto pot = Energy::NonbondedCutg2g<Tspace,Tpairpot>(mcp)
    + Energy::ExternalPressure<Tspace>(mcp)
    + Energy::EquilibriumEnergy<Tspace>(mcp);// + myenergy<Tspace>(mcp);

  auto eqenergy = &pot.second;

  DebyeHuckel dh(mcp);

  //save(pot.first.first.pairpot, atom["ASP"].id, atom["ASP"].id, "pot.dat");
  //save(dh, atom["ASP"].id, atom["ASP"].id, "pot_dh.dat");

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

  Analysis::RadialDistribution<> rdf(0.5);
  Analysis::ChargeMultipole mpol;
  Scatter::DebyeFormula<Scatter::FormFactorUnity<>> debye(mcp);

  spc.load("state");

  Move::Isobaric<Tspace> iso(mcp,pot,spc);
  Move::TranslateRotate<Tspace> gmv(mcp,pot,spc);
  Move::AtomicTranslation<Tspace> mv(mcp,pot,spc);
  Move::SwapMove<Tspace> tit(mcp,pot,spc,*eqenergy);
  if (inPlane)
    for (auto &m : pol)
      gmv.directions[ m.name ]=Point(1,1,0);
  
  gmv.mpi=&mpi;
  tit.mpi=&mpi;

  sys.init( Energy::systemEnergy(spc,pot,spc.p) );    // Store total system energy

  if (mpi.isMaster())
    cout << atom.info() << spc.info() << pot.info() << tit.info()
      << textio::header("MC Simulation Begins!");

  vector<Point> cm_vec; // vector of mass centers

  MCLoop loop(mcp);
  while ( loop.macroCnt() ) {  // Markov chain
    while ( loop.microCnt() ) {
      int k,i=rand() % 2;
      switch (i) {
        case 1:
          k=pol.size();
          while (k-->0) {
            gmv.setGroup( pol[ rand() % pol.size() ] );
            sys+=gmv.move();
          }

          // sample g(r)
          for (auto i=pol.begin(); i!=pol.end()-1; i++)
            for (auto j=i+1; j!=pol.end(); j++)
              rdf( spc.geo.dist(i->cm,j->cm) )++;

          // sample S(q)
          if (slp_global()>1.95) {
            cm_vec.clear();
            for (auto &i : pol)
              cm_vec.push_back(i.cm);
            debye.sample(cm_vec);
          }

          break;
        case 2:
          sys+=iso.move();
          break;
        case 0:
          sys+=tit.move();
          mpol.sample(pol,spc);
          break;
        case 3:
          mv.setGroup(salt);
          sys+=mv.move();
          break;
      }
      if (slp_global()<-0.001 ) {
        xtc.setbox( spc.geo.len );
        xtc.save("traj.xtc", spc);
      }
    } // end of micro loop

    sys.checkDrift( Energy::systemEnergy(spc,pot,spc.p) ); // detect energy drift
    if (mpi.isMaster())
      cout << loop.timing();

  } // end of macro loop

  if (mpi.isMaster()) {
    cout << tit.info() + loop.info() + sys.info() + gmv.info() + mv.info()
    + iso.info() + mpol.info() << endl;
    
    int cnt=0;
    std::ofstream f;
    f.open("debye"+std::to_string(cnt)+".xyz");
    if (f) {
      for (auto m : cm_vec)
        f << "H " << m.transpose() << endl;
      f.close();
    }
    
    rdf.save("rdf_p2p.dat");
    FormatPQR::save("confout.pqr", spc.p, spc.geo.len);
    spc.save("state");
    mcp.save("mdout.mdp");
    debye.save("debye.dat");
  }
}
#endif
