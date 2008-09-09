#include "../container.h"
#include "../point.h"
#include "../slump.h"
#include "../group.h"
#include "../potentials.h"
#include "../species.h"


using namespace std;

int main() {
  cell c(100.);         // Spherical cell with radius 100
  group g1,g2;          // Some groups.
  pot_setup potcfg;
  pot_lB=7.1;           // Specify Bjerrum length
  interaction<pot_coulomb> pot(potcfg); // Specify coulomb potential
  



  //double systemenergy = pot.energy(p);
  //double between2groups = pot.energy(p, g1, g2);
}
