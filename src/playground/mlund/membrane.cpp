#include <faunus/faunus.h>
#include <faunus/membrane.h>

using namespace Faunus;
using namespace Faunus::Potential;

// custom combination of two pair-potentials. T2 if envoked only between certain particles.
template<class T1, class T2>
class CustPairPot : public CombinedPairPotential<T1,T2> {
  private:
    typedef CombinedPairPotential<T1,T2> Tbase;
  public:
    particle::Tid id; // T2 potential added to particle pairs of this type only!
    CustPairPot(InputMap &in) : Tbase(in), id(atom["TL"].id) {}
    inline double operator() (const particle &a, const particle &b, double r2) const FOVERRIDE {
      double u=Tbase::first(a,b,r2);
      if (a.id==id)
        if (b.id==id)
          u+=Tbase::second(a,b,r2);
      return u;
    }
};

typedef Geometry::Cuboid Tgeometry;   // specify geometry - here cube w. periodic boundaries
//typedef CustPairPot<WeeksChandlerAndersen,CosAttract> Tpairpot;
typedef CombinedPairPotential<DebyeHuckel,WeeksChandlerAndersen> Tdhwca;
typedef CustPairPot<Tdhwca,CosAttract> Tpairpot;

int main() {
  cout << textio::splash();           // show faunus banner and credits

  InputMap mcp("membrane.input");     // open user input file
  MCLoop loop(mcp);                   // class for handling mc loops
  FormatPQR pqr;                      // PQR structure file I/O
  FormatAAM aam;                      // AAM structure file I/O
  EnergyDrift sys;                    // class for tracking system energy drifts
  UnitTest test(mcp);                 // class for unit testing

  // Energy functions and space
  Energy::Hamiltonian pot;
  auto nonbonded = pot.create( Energy::Nonbonded<Tpairpot,Tgeometry>(mcp) );
  auto bonded    = pot.create( Energy::Bonded() );
  Space spc( pot.getGeometry() );

  // Markov moves and analysis
  Move::TranslateRotate gmv(mcp,pot,spc);
  Move::AtomicTranslation mv(mcp, pot, spc);
  Move::Pivot pivot(mcp, pot, spc);
  Move::Reptation rep(mcp, pot, spc);
  Analysis::PolymerShape shape;
  Analysis::RadialDistribution<> rdf(0.2);

  // Lipid energyfield
  double sigma = mcp.get<double>("lipid_sigma", 10);
  double epsilon = mcp.get<double>("lipid_epsilon", 1);
  double rho = 1 / mcp.get<double>("lipid_areapercharge",0);
  double headtail_k=0.5*10*epsilon/(sigma*sigma);
  double headtail_req=4*sigma;
  double fene_k=30*epsilon/(sigma*sigma);
  double fene_rmax=1.5*sigma;
  std::shared_ptr<Potential::Harmonic> headtailpot( new Potential::Harmonic(headtail_k,headtail_req) );
  std::shared_ptr<Potential::FENE> fenepot( new Potential::FENE(fene_k, fene_rmax) );

  // Custom Weeks-Chandler-Andersen parameters
  nonbonded->pairpot.first.second.customSigma(atom["HD"].id, atom["HD"].id, 0.95*sigma);
  nonbonded->pairpot.first.second.customSigma(atom["HD"].id, atom["TL"].id, 0.95*sigma);
  nonbonded->pairpot.first.second.customSigma(atom["TL"].id, atom["TL"].id, sigma);
  nonbonded->pairpot.first.second.customEpsilon(atom["HD"].id, atom["HD"].id, epsilon);
  nonbonded->pairpot.first.second.customEpsilon(atom["HD"].id, atom["TL"].id, epsilon);
  nonbonded->pairpot.first.second.customEpsilon(atom["TL"].id, atom["TL"].id, epsilon);

  // Add lipids
  int N=mcp.get("lipid_N",0);
  vector<GroupMolecular> pol(N);   // vector of polymers
  string polyfile = mcp.get<string>("lipid_file", "");
  for (auto &g : pol) {                                  // load polymers
    aam.load(polyfile);
    Geometry::FindSpace f;
    f.find(*spc.geo, spc.p, aam.particles() );           // find empty spot in particle vector
    g = spc.insert( aam.particles() );                   // insert into space
    g.name="lipid";
    spc.enroll(g);
    bonded->add(g.front(), g.front()+1, fenepot);
    bonded->add(g.front()+1, g.front()+2, fenepot);
    bonded->add(g.front(), g.front()+2, headtailpot);

    // dope head groups w. charge
    cout << "# Surface charge density: " << 1/rho << " A^2/charge\n";
    double area=2 * nonbonded->geometry.len.x * nonbonded->geometry.len.y;
    double p = rho*area/N;
    if (slp_global.randOne()<p) {
      spc.p[ g.front() ].charge=-1;
      spc.trial[ g.front() ].charge=-1;
    }
  }

  Group allpol( pol.front().front(), pol.back().back() );// make group w. all polymers

  spc.load("state");                                     // load old config. from disk (if any)
  pqr.save("initial.pqr", spc.p);

  sys.init( Energy::systemEnergy(spc,pot,spc.p)  );      // store initial total system energy

  cout << atom.info() << spc.info() << pot.info() << textio::header("MC Simulation Begins!");

  while ( loop.macroCnt() ) {  // Markov chain 
    while ( loop.microCnt() ) {
      int k,i=slp_global.rand() % 2;
      switch (i) {
        case 0:
          mv.setGroup(allpol);
          sys+=mv.move( allpol.size() ); // translate lipid monomers
          for (auto &g : pol)
            g.setMassCenter(spc);
          break;
        case 1:
          k=pol.size();
          while (k-->0) {
            gmv.setGroup( pol[ slp_global.rand() % pol.size() ] );
            sys+=gmv.move();  // translate/rotate polymers
          }
          break;
        case 2:
          k=pol.size();
          while (k-->0) {
            pivot.setGroup( pol[ slp_global.rand() % pol.size() ] );
            sys+=pivot.move();// pivot move
          }
          break;
      }

      for (auto i=pol.begin(); i!=pol.end()-1; i++)
        for (auto j=i+1; j!=pol.end(); j++)
          rdf( spc.geo->dist(i->cm,j->cm) )++;// polymer mass-center distribution function

    } // end of micro loop

    sys.checkDrift(Energy::systemEnergy(spc,pot,spc.p)); // compare energy sum with current

    // save to disk
    rdf.save("rdf.dat");
    pqr.save("confout.pqr", spc.p);
    spc.save("state");

    cout << loop.timing();

  } // end of macro loop

  // perform unit tests
  gmv.test(test);
  mv.test(test);
  sys.test(test);

  // print information
  cout << loop.info() << sys.info() << mv.info() << gmv.info() << rep.info()
    << pivot.info() << shape.info() << test.info();

  return test.numFailed();
}
