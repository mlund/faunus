#include <faunus/faunus.h>
using namespace Faunus;

typedef Geometry::cuboid Tgeometry;                // select simulation geometry
typedef Potential::coulomb_lj<Tgeometry> Tpairpot; // select particle-particle pairpotential

template<class T>
class distributions {
  private:
    typedef xytable<double,T> Ttable;
    typedef std::map<string,Ttable> Tmap;
    double xmin, xmax, dx;
    Tmap dmap;
  public:
    distributions(double min, double max, double delta) {
    }
    void add(string name, double val) {
      dmap[name]+=val;
    }
};

int main() {
  cout << textio::splash();
  distributions<double> dst(0,100,0.5);
  //dst.add("Utot",0);
  atom.includefile("atomlist.inp");    // load atom properties
  inputfile in("bulk.inp");            // read user input
  mcloop loop(in);                     // class for handling mc loops
  iopqr pqr;                           // PQR structure file I/O
  energydrift sys;                     // class for tracking system energy drifts

  Tgeometry geo(in);                   // simulation geometry
  space spc(geo);                      // generate simulation space
  Energy::nonbonded<Tpairpot> pot(in); // energy calculation class

  // Handle particles
  atom["NA"].dp=20.;                   // Displacement parameter
  atom["CL"].dp=80.;
  spc.insert("NA",100);                // Insert particles into space
  spc.insert("CL",100);
  group salt(0,199);                   // Define salt range
  salt.name="Salt particles";
  Move::translate_particle mv(in, pot, spc);  // Particle move class
  mv.igroup = spc.enroll(salt);               // Enroll salt in space and select it for particle moves

  spc.load("space.state");
  sys.init( pot.all2all(spc.p) );

  cout << atom.info() << spc.info() << pot.info() << mv.info() << pot.pair.info() << endl;

  cout << textio::header("MC Simulation Begins!");

  while ( loop.macroCnt() ) {  // Markov chain 
    while ( loop.microCnt() ) {
      sys+=mv.move();
    }
    sys.checkdrift( pot.all2all(spc.p) );
    cout << loop.timing();
  }

  pqr.save("confout.pqr", spc.p);
  spc.save("space.state");
  cout << mv.info() << sys.info() << loop.info();
}
