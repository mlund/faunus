#include <iostream>
#include "../reactionfield.h"
#include "../point.C"
#include "../interact.C"
#include "../slump.C"
#include "../group.C"

using namespace std;

int main() {
  double a=15.;
  interact_rf rf(7.14, a, 80, 2);

  particle p1,p2;

  p1.x=10;
  p1.charge=+1.0;
  p2.z=10;

  //cout << "rfield  = " << rf.lB*rf.phi(p1,p2,10) << endl;
  //cout << "coulomb = " << rf.lB*p1.charge / p1.dist(p2) << endl;

  //calculate interaction between two charges at different
  //angles around a sphere
  double pi=acos(-1.);
  spherical s1(52.,0,0), s2(52.,pi,0);
  p1.clear();
  p2.clear();
  p1.charge=+1;
  p2.charge=+1;
  p1=s1.cartesian();
  p2=s2.cartesian();
  for (s2.theta=0.1; s2.theta<pi; s2.theta+=0.01) {
    interact_rf pot(7.14, a, 80, 2);
    p2=s2.cartesian();
    double u=pot.lB*pot.eo * ( pot.energy_rf(p1,p2) );
    //cout << s2.theta << " " << u << " "
    //     << pot.lB*pot.energy(p1,p2) << endl;
  };
  //return 0;


  //Calculate self-energy of charge as a function of
  //the distance to the low dielectric cavity
  p1.clear();
  p1.charge=1;
  p1.radius=2.;
  for (p1.z=a+p1.radius; p1.z<30.; p1.z+=0.01)
    cout << p1.z << " " << rf.lB*rf.eo*rf.selfenergy(p1) << endl;
};
