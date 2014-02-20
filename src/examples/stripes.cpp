/*
For more info: doi:10.1038/nmat820
*/

#include <faunus/faunus.h>
using namespace Faunus;                   // use Faunus namespace
typedef Space<Geometry::Cuboid> Tspace;   // Type of simulation space

namespace Faunus {
  namespace Potential {
    // let's make a new pair potential
    struct CoreShell : public PairPotentialBase {
      double a1,a2,eps;
      inline CoreShell(InputMap& in) {
        name="Coreshell";
        a1 = pow( 2*in("core_radius",1.0), 2);
        a2 = pow( 2*in("shell_radius",2.5), 2);
        eps = in("epsilon", 0.2);
      }
      template<class Tparticle>
        double operator() (const Tparticle &a, const Tparticle &b, double r2) {
          if (r2>a2) return 0;
          if (r2<a1) return pc::infty;
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

  for (auto &i : spc.p)                   // Place particles
    i.z()=0;                              // in XY plane
  spc.trial=spc.p;

  spc.load("state");

  Move::AtomicTranslation<Tspace> mv(in,pot,spc);// particle move class
  mv.dir=Point(1,1,0);                    // move only in xy plane
  mv.setGroup(salt);                      // move class acts on salt group
  mv.move( in("steps",1e3) );             // move salt randomly 'steps' times

  cout << spc.info() + pot.info() + mv.info(); // final information

  spc.save("state");
  FormatPQR::save("stripes.pqr", spc.p);
}
