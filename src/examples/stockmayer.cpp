#include <faunus/faunus.h>
using namespace Faunus;                               // use Faunus namespace
using namespace Faunus::Move;

typedef Geometry::Cuboid Tgeo;                        // select simulation geometry and pair potential
typedef Potential::LennardJones Tpair;

int main() {
  ::atom.includefile("stockmayer.json");              // load atom properties
  InputMap in("stockmayer.input");                    // open parameter file for user input
  Energy::Nonbonded<Tpair,Tgeo> pot(in);              // create Hamiltonian, non-bonded only
  Space spc( pot.getGeometry() );                     // create simulation space, particles etc.
  GroupAtomic sol(spc, in);                           // group for particles
  MCLoop loop(in);                                    // class for mc loop counting
  Analysis::RadialDistribution<> rdf(0.2);            // particle-particle g(r)

  Move::AtomicTranslation trans(in, pot, spc);        // particle move class
  Move::AtomicRotation rot(in, pot, spc);             // particle move class
  //PolarizeMove<AtomicTranslation> trans(in,pot,spc);
  //PolarizeMove<AtomicRotation> rot(in,pot,spc);
  trans.setGroup(sol);                                // tells move class to act on sol group
  rot.setGroup(sol);                                  // tells move class to act on sol group

  while ( loop.macroCnt() ) {                         // Markov chain 
    while ( loop.microCnt() ) {
      if (slp_global() > 0.5)
        trans.move( sol.size() );                     // translate
      else
        rot.move( sol.size() );                       // rotate
    }
  }

  rdf.save("gofr.dat");                               // save g(r) to disk
  std::cout << spc.info() + pot.info() + trans.info() + rot.info(); // final info
}
