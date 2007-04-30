#include "container.h"

//----------- CONTAINER -----------------
group container::insert(particle a, short n) {
  hardsphere hs;
  group g;
  if (n==0)
    return g;
  g.beg=p.size();
  for (int i=0; i<n; i++) {
    g.end=push_back(a)-1;
    while (hs.overlap(p, g.end)==true) {
      randompos(p[g.end]);
      trial[g.end] = p[g.end];
    }
  }
  return g;
}


//----------- CELL ----------------------
cell::cell(float radius) {
  r = radius; 
  r2 = r*r; 
  diameter = 2*r; 
  volume = (4./3.)*acos(-1)*r*r*r;
}
void cell::randompos(point &p) {
  double l=r2+1;
  while (l>r2) {
    p.x = slp.random_half()*diameter;
    p.y = slp.random_half()*diameter;
    p.z = slp.random_half()*diameter;
    l=p.x*p.x+p.y*p.y+p.z*p.z;
  }
}
//----------- BOX --------------------------
box::box(float sidelength) {
  len = sidelength;
  volume = len*len*len;
}

void box::randompos(point &p) {
}
