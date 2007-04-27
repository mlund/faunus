#include "container.h"

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
