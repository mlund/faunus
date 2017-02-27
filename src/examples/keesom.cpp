#include <faunus/faunus.h>

using namespace Faunus;

// Keesom potential w. a Bjerrum length of 7 angstroms
double keesom(double mu, double r) {
  return -pow(7.*mu*mu,2) / (3*pow(r,6)); // in kT
}

int main() {
  DipoleParticle a;            // create a dipole particle
  a.clear();                   // ...and make sure it is empty 
  a.setMu(Point(1,0,0));              // initial dipole direction (unit vector)
  a.setMuscalar(2.2_Debye);      // ...and  magnitude (see also "_eA", "_Cm")
  auto b=a;                    // duplicate particle...
  b.x() = 0.9_angstrom;        // ...and separate by 0.9 angstrom

  Potential::DipoleDipole u(298.15_K, 80.); // pair potential
  RandomTwister<> ran;         // random number generator

  while (b.x() < 0.4_nm) {     // loop over separations
    Average<double> Q;         // partition function (just an average)
    auto r = a-b;              // separation vector
    for (int i=0; i<1e5; i++) {// loop over angular space...
      a.mu().ranunit( ran );     // ...via random dipole directions
      b.mu().ranunit( ran );     //    - / / -
      Q += exp( -u(a,b,r) );   // Boltzmann weighted dipole-dipole energy
    }
    printf("%9.2f %9.5f %9.5f\n",
        r.norm(), -log(Q), keesom(a.muscalar(), r.norm()) );

    b.x() += 0.1_angstrom;     // increase separation
  }
}

/**  @page example_keesom Example: Angularly Averaged Dipole-Dipole Energy

 This example demonstrates how to construct particles
 using arbitrary units for i.e. dipole moment and
 distances. For more info, see namespace `ChemistryUnits`.
 It also shows how to calculate the inter-particle energy
 which, by brute force, is used to evaluate the angularly averaged
 potential of mean force (interaction free energy) between a pair
 of dipoles,
 @f[
   \beta w(r)_{exact}
   = -\ln \langle
   \exp{(-\beta u_{dd}(\boldsymbol{\mu}_1,\boldsymbol{\mu}_2 )} \rangle
 @f]
 where the square brackets denote an average over angular space.
 This *exact* numerical result will asymptotically approach
 the approximate Keesom potential at long separations,
 @f[
   \beta w(r)_{keesom} = - \frac{ (\lambda_B \mu_1\mu_2)^2 }{3r^6}
 @f]
 where @f$ \lambda_B @f$ is the Bjerrum length.
 
 @includelineno examples/keesom.cpp
*/
