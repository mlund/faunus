#include <faunus/faunus.h>

using namespace Faunus;

typedef Geometry::Cuboid Tgeometry;                // select simulation geometry
typedef Potential::CoulombSR<Tgeometry, Potential::Coulomb, Potential::LennardJones> Tpairpot;
//typedef Potential::CoulombSR<Tgeometry, Potential::Coulomb, Potential::HardSphere> Tpairpot;

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
  FormatPQR pqr;                           // PQR structure file I/O
  EnergyDrift sys;                     // class for tracking system energy drifts

  Energy::Hamiltonian pot;
  auto nonbonded = pot.create( Energy::Nonbonded<Tpairpot>(mcp) );
  auto bonded    = pot.create( Energy::Bonded() );
  Space spc( pot.getGeometry() );

  // Handle particles
  GroupAtomic salt(spc, mcp);
  salt.name="Salt";
  Move::ParticleTranslation mv(mcp, pot, spc);  // Particle move class
  mv.setGroup(salt);
  if (salt.id==Group::ATOMIC)
    cout << "ATOMIC!!\n";

  GroupMolecular test;
  test.beg=spc.p.size();
  spc.insert("Na", 2);
  test.end=spc.p.size()-1;
  test.name="Test";
  spc.p[test.beg].x=0;
  spc.p[test.beg].y=0;
  spc.p[test.beg].z=0;
  spc.p[test.end].x=2;
  spc.p[test.end].y=2;
  spc.p[test.end].z=2;
  spc.trial=spc.p;
  spc.enroll(test);

  spc.load("space.state");
  salt.setMassCenter(spc);
  test.setMassCenter(spc);

  //bonded->bonds.add(0,1, Potential::Harmonic(0.2, 10.0));
  //bonded->bonds.add(1,2, Potential::Harmonic(0.3,  5.0));

  // Particle titration
  Move::SwapMove tit(mcp,pot,spc);
  Move::Isobaric iso(mcp,pot,spc);
  Move::RotateGroup gmv(mcp,pot,spc);

  // Widom particle insertion
  Analysis::Widom widom(spc, pot);
  widom.addGhost(spc);

  //FastaSequence fasta;
  //Group protein1 = fasta.insert( "AAK", spc, bonded->bonds );
  
  #define UTOTAL \
      pot.g_external(spc.p, test)\
    + pot.g_internal(spc.p, salt)  + pot.g_external(spc.p, salt)\
    + pot.g2g(spc.p,salt,test) + pot.external()

  //#define UTOTAL pot.g_internal(spc.p, salt)  + pot.g_external(spc.p, salt) + pot.external()

  sys.init( UTOTAL );

  cout << atom.info() << spc.info() << pot.info() << mv.info()
       << textio::header("MC Simulation Begins!");

  while ( loop.macroCnt() ) {  // Markov chain 
    while ( loop.microCnt() ) {
      sys+=mv.move();
      //if (slp_global.randOne()>0.9)
      //  widom.sample();
      sys+=tit.move();
      gmv.setGroup(test);
      sys+=gmv.move();
      sys+=iso.move();
    }
    sys.checkDrift( UTOTAL );
    cout << loop.timing();
  }

  pqr.save("confout.pqr", spc.p);
  spc.save("space.state");

  cout << mv.info() << sys.info() << loop.info() << widom.info()
       << tit.info() << gmv.info() << iso.info();
}
