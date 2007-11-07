#include "../potentials.h"
#include <iostream>

using namespace std;

int main() {
  interaction p(7.1,1.);

  particle a,b;
  a.z=10;
  b.z=-10;
  a.charge=1;
  b.charge=-1;
  a.radius=1.5;
  b.radius=1.5;

  for (double r=0; r<100; r+=0.25) {
    a.z=r;
    b.z=-r;
    cout << a.z-b.z << " " << p.f*p.pairpot(a,b) << endl;
  }
};
