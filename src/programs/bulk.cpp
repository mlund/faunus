#include <faunus/faunus.h>
#include <faunus/widom.h>

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
  //dst.add("Utot",0);
  atom.includefile("atomlist.inp");    // load atom properties
  InputMap mcp("bulk.inp");
  MCLoop loop(mcp);                    // class for handling mc loops
  iopqr pqr;                           // PQR structure file I/O
  energydrift sys;                     // class for tracking system energy drifts

  Energy::Nonbonded<Tpairpot> pot(mcp); // energy calculation class
  Space spc( pot.getGeometry() );

  // Handle particles
  atom["NA"].dp=20.;                   // Displacement parameter
  atom["CL"].dp=80.;
  spc.insert("NA",100);                // Insert particles into space
  spc.insert("CL",100);
  group salt(0,199);                   // Define salt range
  salt.name="Salt particles";
  spc.enroll(salt);
  Move::ParticleTranslation mv(mcp, pot, spc);  // Particle move class
  mv.setGroup(salt);

  spc.load("space.state");
  sys.init( pot.all2all(spc.p) );

  Widom wid(spc, pot);
  wid.addGhost(spc);


  cout << atom.info() << spc.info() << pot.info() << mv.info() << endl;

  cout << textio::header("MC Simulation Begins!");

  while ( loop.macroCnt() ) {  // Markov chain 
    while ( loop.microCnt() ) {
      sys+=mv.move();
      wid.sample();
    }
    sys.checkdrift( pot.all2all(spc.p) );
    cout << loop.timing();
  }

  pqr.save("confout.pqr", spc.p);
  spc.save("space.state");
  cout << mv.info() << sys.info() << loop.info() << wid.info();
}
