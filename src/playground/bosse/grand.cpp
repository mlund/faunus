#include <faunus/faunus.h>
#include <faunus/ewald.h>
using namespace Faunus;
using namespace Faunus::Potential;

typedef Space<Geometry::Cuboid,PointParticle> Tspace;

// Define custom pair potential
class Cardinaux : public PairPotentialBase {
  private:
    string _brief() {
      return name+": a=" + std::to_string(alpha);
    }
    int alpha,alphahalf;
    PairMatrix<double> eps;
  public:
    Cardinaux(InputMap &in, string pfx="cardinaux") {
      name="Cardinaux";
      alpha=in.get<int>(pfx+"_alpha", 90);
      alphahalf=alpha/2;
      for (auto i : atom.list)
        for (auto j : atom.list)
          eps.set(i.id, j.id, 4*pc::kJ2kT(
                sqrt(atom[i.id].eps*atom[j.id].eps) // epsilon mixing
                ));
    }
    template<class Tparticle>
      double operator() (const Tparticle &a, const Tparticle &b, double r2) const {
        double s=a.radius+b.radius;
        s=s*s/r2; // ^2
        s=__builtin_powi(s,alphahalf); // = (sigma/r)^a [gcc/clang specific pow()!]
        return eps(a.id,b.id)*(s*s - s);
      }
    template<typename Tparticle>
      Point force(const Tparticle &a, const Tparticle &b, double r2, const Point &p) {
        double s=a.radius+b.radius;
        s=s*s/r2; //^2
        s=__builtin_powi(s,alphahalf); // = (s/r)^a
        return alpha*eps(a.id,b.id)*s*(2*s-1)/r2 * p; // extra division can be avoided!
      }
    string info(char w) { return "  "+_brief()+"\n"; }
};

typedef CombinedPairPotential<Coulomb,Cardinaux> Tpairpot;

int main() {
  InputMap mcp("grand.input");                  // open user input file
  Tspace spc(mcp);                              // simulation space
  Energy::Nonbonded<Tspace,Tpairpot> pot(mcp);  // hamiltonian
  pot.setSpace(spc);                            // share space w. hamiltonian

  Group salt;                                   // group for salt particles
  salt.addParticles(spc, mcp);
  spc.load("state",Tspace::RESIZE);             // load old config. from disk (if any)

  Analysis::VirialPressure virial;
  Move::GrandCanonicalSalt<Tspace> gc(mcp,pot,spc,salt);
  Move::AtomicTranslation<Tspace> mv(mcp,pot,spc);
  mv.setGroup(salt);

  EnergyDrift sys;                              // class for tracking system energy drifts
  sys.init(Energy::systemEnergy(spc,pot,spc.p));// store initial total system energy

  cout << atom.info() + spc.info() + pot.info() + "\n";

  save(pot.pairpot, atom["Na"].id, atom["Cl"].id, "na-cl.dat");
  save(pot.pairpot, atom["Na"].id, atom["M"].id, "na-m.dat");
  save(pot.pairpot, atom["M"].id, atom["M"].id, "m-m.dat");

  MCLoop loop(mcp);                             // class for handling mc loops
  while ( loop.macroCnt() ) {
    while ( loop.microCnt() ) {
      if (slp_global() < 0.5)
        sys+=mv.move( salt.size() );            // translate salt
      else 
        sys+=gc.move( salt.size()/2 );          // grand canonical exchange
      virial.sample(spc, pot);                  // sample virial
    }                                           // end of micro loop
    sys.checkDrift(Energy::systemEnergy(spc,pot,spc.p)); // calc. energy drift
    cout << loop.timing();
  }                                             // end of macro loop

  FormatPQR::save("confout.pqr", spc.p);        // PQR snapshot for VMD etc.
  spc.save("state");                            // final simulation state

  cout << loop.info() + sys.info() + mv.info() + gc.info() + virial.info();
}
