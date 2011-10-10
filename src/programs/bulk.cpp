#include <faunus/faunus.h>

using namespace Faunus;

typedef Geometry::Cuboid Tgeometry;                // select simulation geometry
typedef Potential::coulomb_lj<Tgeometry> Tpairpot; // select particle-particle pairpotential

template<class T>
class distributions {
  private:
    typedef xytable<double,T> Ttable;
    typedef std::map<string,Ttable> Tmap;
    double xmin, xmax, dx;
    Tmap map;
  public:
    distributions(double min, double max, double delta) {
    }
    void add(string name, double val) {
      map[name]+=val;
    }
};

int main() {
  cout << textio::splash();
  distributions<double> dst(0,100,0.5);
  atom.includefile("atomlist.inp");    // load atom properties
  InputMap mcp("bulk.inp");
  MCLoop loop(mcp);                    // class for handling mc loops
  iopqr pqr;                           // PQR structure file I/O
  energydrift sys;                     // class for tracking system energy drifts

  Energy::Nonbonded<Tpairpot> nb(mcp); // non-bonded energy
  Energy::Hamiltonian pot;
  pot.add(nb);

  Space spc( pot.getGeometry() );

  // Handle particles
  Atomic salt(spc, mcp);
  salt.name="Salt particles";
  Move::ParticleTranslation mv(mcp, pot, spc);  // Particle move class
  mv.setGroup(salt);
  spc.load("space.state");

  // Particle titration
  Move::SwapMove tit(mcp,pot,spc);

  // Widom particle insertion
  Analysis::Widom widom(spc, pot);
  widom.addGhost(spc);

  sys.init( mv.totalEnergy() + tit.totalEnergy() );

  cout << atom.info() << spc.info() << pot.info() << mv.info()
       << textio::header("MC Simulation Begins!");

  while ( loop.macroCnt() ) {  // Markov chain 
    while ( loop.microCnt() ) {
      sys+=mv.move();
      if (slp_global.random_one()>0.9)
        widom.sample();
      sys+=tit.move();
    }
    sys.checkdrift( mv.totalEnergy() + tit.totalEnergy()  );
    cout << loop.timing();
  }

  pqr.save("confout.pqr", spc.p);
  spc.save("space.state");
  cout << mv.info() << sys.info() << loop.info() << widom.info();
  cout << tit.info();
}
