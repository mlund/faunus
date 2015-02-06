#include <faunus/faunus.h>

using namespace Faunus;                   // use Faunus namespace
typedef Space<Geometry::Cuboid> Tspace;   // Type of simulation space

namespace Faunus {
  namespace Potential {
    // let's make a new pair potential
    struct CoreShell : public PairPotentialBase {
      double a1,a2,eps;
      inline CoreShell(InputMap& in, string dir="") : PairPotentialBase(dir) {
        name = "Core-shell potential";
        in.cd ( jsondir + "/coreshell" );
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
  InputMap in("stripes.json");            // open parameter file for user input
  Tspace spc(in);                         // simulation space, particles etc.
  Energy::Nonbonded<Tspace,Potential::CoreShell> pot(in);// Hamiltonian, non-bonded only

  spc.load("state");                      // load old configuration if any

  Move::AtomicTranslation<Tspace> mv(in,pot,spc);// particle move class

  for (int i=0; i<1e3; i++)
    mv.move();

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
 potential that takes user input from a json file.

 Run this example from the `examples` directory and for the first run,
 just ignore any potential warnings about infinite energy which is
 caused by initial particle overlap.

 ~~~~~~~~~~~~~~~~~~~
 $ make example_stripes
 $ cd src/examples
 $ ./stripes.py
 ~~~~~~~~~~~~~~~~~~~

 ![Stripes? T=0.18, rho=0.291, N=1000](stripes.jpg)

 stripes.cpp
 =============
 @includelineno examples/stripes.cpp

 stripes.json
 ============
 @includelineno examples/stripes.json

*/
