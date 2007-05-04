#include "container.h"

//----------- CONTAINER -----------------
string container::info() {
  float z=charge();
  ostringstream o;
  o << "# Container:" << endl
    << "#   Number of particles  = " << p.size() << endl
    << "#   Volume (AA^3)        = " << volume << endl
    << "#   Electroneutrality    = " 
    << ((z!=0) ? "NO!" : "Yes") << " "  << z << endl;;
  return o.str();
}
group container::insert(particle::type id, short n) {
  group g;
  hardsphere hs;
  particle a=get(id);
  g.beg=p.size();
  for (unsigned short i=0; i<n; i++) {
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
string cell::info() {
  ostringstream o;
  o << container::info() 
    << "#   Shape                = Spherical" << endl
    << "#   Radius               = " << r << endl;
  return o.str();
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
