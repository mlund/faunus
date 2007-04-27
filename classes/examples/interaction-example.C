#include "../point.C"
#include "../slump.C"
#include "../group.C"
#include "../potentials.C"
#include "../species.C"

using namespace std;

int main() {
  group g1,g2;
  vector<particle> p(1);
  pot_setup potcfg;
  potcfg.lB=7.1;
  potcfg.eps=1.0;
  interaction<pot_coulomb> pot(potcfg);

  double systemenergy = pot.energy(p);
  double between2groups = pot.energy(p, g1, g2);
}
