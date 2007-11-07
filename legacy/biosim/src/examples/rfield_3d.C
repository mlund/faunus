#include <iostream>
#include "../reactionfield.h"
#include "../point.C"
#include "../interact.C"
#include "../slump.C"
#include "../group.C"

using namespace std;

int main() {
  double a=13.01, box=40, phi;
  double pi=acos(-1.);
  interact_rf pot(7.14, a, 80, 20);
  particle p1,p2,p3;
  spherical s1(15.01,0,0), s2(15.01,pi/2.,0);
  p1.charge=+1;
  p2.charge=+1;
  p1 = s1.cartesian();
  p2 = s2.cartesian();

  for (double x=-box; x<box; x+=1.) {
    //cout << endl << x << endl;
    for (double z=-box; z<box; z+=1.) {
      p3.x=x;
      p3.z=z;
      phi = pot.lB*pot.eo * ( pot.phi(p1,p3) + pot.phi(p2,p3) );
      phi = phi / (pot.lB*p1.charge/p1.dist(p3) + pot.lB*p2.charge/p2.dist(p3));
      if (sqrt(x*x+z*z)>0.01) //origo collision
	cout << x << " " << z << " " << phi << endl;
    };
    cout << endl;
  };
  return 0;
};
