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
        a1 = pow( in("core_radius",1.0), 2);
        a2 = pow( in("shell_radius",2.5), 2);
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
  Tspace spc(in);                         // simulation space, particles etc.
  Group salt;                             // group for salt particles
  salt.addParticles(spc,in);              // add according to user input

  for (auto &i : spc.p)                   // place particles
    i.z()=0;                              // in XY plane
  spc.trial=spc.p;                        // trial coords must be in sync

  spc.load("state");                      // load old configuration if any

  Move::AtomicTranslation<Tspace> mv(in,pot,spc);// particle move class
  mv.dir=Point(1,1,0);                    // move only in xy plane
  mv.setGroup(salt);                      // move class acts on salt group
  mv.move( in("steps",1e3) );             // move salt randomly 'steps' times

  spc.save("state");                      // save final state
  FormatPQR::save("stripes.pqr", spc.p);  // save PQR file for i.e. VMD

  cout << spc.info() + pot.info() + mv.info(); // final information
}

/** @page example_stripes Example: Stripe Phase from Isotropic Repulsion
 
 This will simulate core-shell particles on a two-dimensional
 surface. The pair interaction is similar to a hard core
 potential but with a step, mimicking an outer core.
 As has been shown (http://dx.doi.org/10.1038/nmat820),
 this simple system shows remarkable phases with strong ordering.
 This particular example also shows how to create a simple, custom pair
 potential that takes user input from the main input file.

 Run this example from the `examples` directory:

 ~~~~~~~~~~~~~~~~~~~
 $ make examples_stripes
 $ cd src/examples
 $ ./stripes.run
 ~~~~~~~~~~~~~~~~~~~

 ![Stripes? T=0.18, rho=0.291, N=1000](stripes.jpg)

 stripes.cpp
 =============
 @includelineno examples/stripes.cpp

 stripes.run
 =============
 @includelineno examples/stripes.run

*/
