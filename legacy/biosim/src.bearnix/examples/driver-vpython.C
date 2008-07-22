#include <iostream>
#include <vector>
#include "vpython.C"
#include "group.C"
#include "slump.C"
#include "point.C"

using namespace std;

int main() {
  vpython vp;
  slump s;
  group g;

  g.beg = 0;
  g.end = 9;
  vector<particle::particle> p(10);
  
  for (int i=0; i<p.size(); i++) {
    p[i].x=10*s.random_one();
    p[i].y=10*s.random_one();
    p[i].z=10*s.random_one();
    p[i].radius=0.001;
  };
  
  vp.add( p, g );
  vp.save("sletmig.py");
  
};
