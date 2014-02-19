#include <faunus/faunus.h>
using namespace Faunus;                   // use Faunus namespace
typedef Space<Geometry::Cuboid> Tspace;   // Type of simulation space

namespace Faunus {
  namespace Potential {
    struct CoreShell : public PairPotentialBase {
      double a1,a2,eps;
      inline CoreShell(InputMap& in) { name="Hardsphere"; }
      inline string info(char w) { return name; }
      template<class Tparticle>
        double operator() (const Tparticle &a, const Tparticle &b, double r2) {
          if (r2>a2*a2) return 0;
          if (r2<a1*a1) return pc::infty;
          return eps;
        }
    };
  }
}

int main() {
  ::atom.includefile("stripes.json");     // load atom properties
  InputMap in("stripes.input");           // open parameter file for user input
  Energy::Nonbonded<Tspace,Potential::CoreShell> pot(in);// Hamiltonian, non-bonded only
  Tspace spc(in);                         // Simulation space, particles etc.
  Group salt;                             // Group for salt particles
  salt.addParticles(spc,in);              // Add according to user input

  for (auto &i : spc.p)
    i.z()=0;
  spc.trial=spc.p;

  Move::AtomicTranslation<Tspace> mv(in,pot,spc);// particle move class
  mv.dir=Point(1,1,0);
  mv.setGroup(salt);                      // move class acts on salt group
  mv.move(1e5);                           // move salt randomly 100000 times

  cout << spc.info() + pot.info() + mv.info(); // final information

  FormatPQR::save("stripes.pqr", spc.p);
}
