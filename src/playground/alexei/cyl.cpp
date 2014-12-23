#include <faunus/faunus.h>

using namespace Faunus;
using namespace Faunus::Potential;

double q2q(double chargeA, double chargeB, double r) {
  return chargeA*chargeB/r;
}

template<bool useIonIon=false, bool useIonDipole=false, bool useDipoleDipole=false, bool useIonQuadrupole=false>
class Multipole : public PairPotentialBase {
  private:
    string _brief() {
      std::ostringstream o;
      o << "Multipole , lB=" << _lB << textio::_angstrom;
      return o.str();          
    }
  protected:
    double _lB;
  public:
    Multipole(InputMap &in) {
      name="Multipole";
      pc::setT ( in.get<double>("temperature", 298.15,
            "Absolute temperature (K)") );
      double epsilon_r = in.get<double>("epsilon_r",80.,
          "Dielectric constant");
      _lB = pc::lB(epsilon_r);
    }
    template<class Tparticle>
      double operator()(const Tparticle &a, const Tparticle &b, const Point &r) {
        double U_total = 0;
        if(useIonIon == true) U_total += q2q(a.charge,b.charge,r.norm());
        if(useIonDipole == true) U_total += q2mu(a.charge*b.muscalar,b.mu,b.charge*a.muscalar,a.mu,r);
        if(useDipoleDipole == true) U_total += mu2mu(a.mu,b.mu, a.muscalar*b.muscalar, r);
        if(useIonQuadrupole == true) U_total += q2quad(a.charge, b.theta,b.charge, a.theta,r);
        return _lB*U_total;
      }

    template<class Tparticle>
      Point field(const Tparticle &p, const Point &r) {
        return Point();
      }

    template<class Tparticle>
      double fieldEnergy(const Tparticle &p, const Point &E) {
        return 0;
      }

    string info(char w) {
      using namespace textio;
      std::ostringstream o;
      o << pad(SUB,w,"Temperature") << pc::T() << " K" << endl
        << pad(SUB,w,"Bjerrum length") << _lB << " "+angstrom << endl;
      return o.str();
    }
};

// class for calculating the system energy
template<class Tspace, class Tenergy, class Tpvec>
double mySystemEnergy(Tspace &spc, Tenergy &pot, const Tpvec &p) {
  pot.setSpace(spc); // ensure pot geometry is in sync with spc
  double u = pot.external(p);
  for (auto g : spc.groupList())
    if (!g->isMolecular())
      u += pot.g_internal(p, *g);
  for (size_t i=0; i<spc.groupList().size()-1; i++)
    for (size_t j=i+1; j<spc.groupList().size(); j++)
      u += pot.g2g(p, *spc.groupList()[i], *spc.groupList()[j]);
  return u;
}

typedef Space<Geometry::PeriodicCylinder, DipoleParticle> Tspace;
typedef CombinedPairPotential<HardSphere, Multipole<true,true,true> > Tpairpot;

int main(int argc, char** argv) {
  InputMap mcp("cyl.input");
  MCLoop loop(mcp);                    // class for handling mc loops
  FormatXTC xtc(1000);                 // XTC gromacs trajectory format
  EnergyDrift sys;                     // class for tracking system energy drifts
  UnitTest test(mcp);

  bool movie = mcp.get("movie", true);

  Tspace spc(mcp);
  auto pot = Energy::NonbondedVector<Tspace,Tpairpot>(mcp);
  pot.setSpace(spc);

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
    FormatPQR::load(file,v);

    Geometry::FindSpace fs;
    if (inPlane)
      fs.dir=Point(0,0,1);
    fs.find(spc.geo,spc.p,v);
    pol[i] = spc.insert(v);
    pol[i].name=file+std::to_string(i);
    spc.enroll( pol[i] );
  }
  Group allpol( pol.front().front(), pol.back().back() );

  // manually set dipole
  double muscalar = mcp.get("muscalar", 1.0) * 1.0_Debye;
  spc.p[1].mu = {1,0,0};
  spc.p[3].mu = {1,0,0};
  spc.p[1].muscalar = muscalar;
  spc.p[3].muscalar = muscalar;
  spc.trial = spc.p;

  // Add atomic species
  Group salt;
  salt.addParticles(spc, mcp);
  salt.name="Atomic Species";
  spc.enroll(salt);

  Analysis::LineDistribution<> rdf(0.5);

  spc.load("state");

  Move::AtomicTranslation<Tspace> mv(mcp, pot, spc);
  Move::TranslateRotate<Tspace> gmv(mcp,pot,spc);
  if (inPlane)
    for (auto &m : pol)                                                       
      gmv.directions[ m.name ]=Point(0,0,1);    

  sys.init( mySystemEnergy(spc,pot,spc.p) );    // Store total system energy

  cout << atom.info() << spc.info() << pot.info() 
    << textio::header("MC Simulation Begins!");

  cout << "# Dipole moment 1: " << spc.p[1].muscalar << " eA (" << spc.p[1].mu.transpose() << ")\n";
  cout << "# Dipole moment 2: " << spc.p[3].muscalar << " eA (" << spc.p[3].mu.transpose() << ")\n";

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
          mv.setGroup(salt);
          sys+=mv.move();
          break;
      }

      double xi = slump(); // random number [0,1)

      if ( movie==true && xi>0.995 ) {
        xtc.setbox( 1000. );
        xtc.save("traj.xtc", spc.p);
      }
    } // end of micro loop

    sys.checkDrift( mySystemEnergy(spc,pot,spc.p) ); // detect energy drift
    cout << loop.timing();

    spc.save("state");
    mcp.save("mdout.mdp");
    rdf.save("rdf_p2p.dat");
    FormatPQR::save("confout.pqr", spc.p);

  } // end of macro loop

  cout << loop.info() + sys.info() + gmv.info() + mv.info();

}
