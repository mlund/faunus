#include <faunus/faunus.h>

using namespace Faunus;

int main() {
  DipoleParticle a;                   // create a dipole particle
  a.clear();                        
  a.mu = {1, 0, 0};                   // dipole moment unit vector
  a.muscalar = 2.2_Debye;             // set dipole moment scalar
  auto b=a;                           // duplicate particle...
  b.x() = 1.0_angstrom;               // ...and separate by 1 angstrom

  Potential::DipoleDipole pot(298.15_K, 80.); // pair potential
  RandomTwister ran;                  // random number generator

  while (b.x() < 1.0_nm) {            // loop over separations
    Average<double> u;                // angularly averaged energy
    double r2 = (a-b).squaredNorm();  // squared distance
    for (int i=0; i<1e6; i++) {
      a.mu.ranunit( ran );            // random dipole direction
      b.mu.ranunit( ran );            //    - / / -
      u += pot(a,b,r2);               // calc. dipole-dipole energy and add to avg.
    }
    std::cout << sqrt(r2) << " " << u << "\n"; // angstrom and kT
    b.x() += 0.05_nm;               // increase separation
  }
}

/**  @page example_keesom Example: Angularly Averaged Dipole-Dipole Energy

 This example demonstrates how to construct particles
 using arbitrary units for i.e. dipole moment and
 distances. For more info, see namespace `ChemistryUnits`.
 It also shows how to calculate inter-particle energies
 which will be used to evaluate the angularly averaged
 potential of mean force between a pair of dipoles.
 The exact numerical result will asymptotically approach
 the approximate Keesom potential.
 
 @includelineno examples/keesom.cpp

*/
