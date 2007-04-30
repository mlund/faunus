#include <iostream>
#include "container.C"
#include "point.C"
#include "hardsphere.C"
#include "species.C"
#include "space.C"
#include "slump.C"
#include "group.C"
#include "potentials.C"
#include "markovmove.C"

using namespace std;

typedef pot_coulomb T_pairpot;

int main() {
  space s;
  cell cell(100.);
  canonical canon;
  pot_setup cfg;
  interaction<pot_coulomb> pot(cfg);

  saltmove sm(canon, s, pot, cell);

};

