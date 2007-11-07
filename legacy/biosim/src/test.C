#include <iostream>
#include "reactionfield.h"
#include "physconst.h"
#include "point.h"

using namespace std;

int main() {

  double a=10.;
  double eo=80,ei=80;
  
  physconst phys;
  interact_rf pot(phys.beta_ecf, a, eo, ei);
  particle p1,p2;

  p1.charge=p2.charge=-1.0;

  p1.x=6;
  p2.x=25.;

  cout << "phi1 = " << pot.phi(p1,p2,false) << endl;
  cout << "self = " << pot.selfenergy(p2) << endl;
};
