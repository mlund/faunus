#include "point.h"

/*!
 *
 * JUST SOMETHING FROM THE C++ FILE
 *
 */

//Constructor
point::point() {
   clear();
};

void point::clear() { x=y=z=0; };

//! length of vector
double point::len() {
  double l2=x*x+y*y+z*z;
  if (l2!=0)
    return sqrt(l2);
  return 0;
};

void point::operator+=(point p) {
  x += p.x;
  y += p.y;
  z += p.z;
};

point point::operator-() {
  point o;
  o.x=-x; o.y=-y; o.z=-z;
  return o;
};

point point::operator*(point p) {
  point o;
  o.x = p.x * x;
  o.y = p.y * y;
  o.z = p.z * z;
  return o;
};
point point::operator*(double s) {
  point o;
  o.x = x*s;
  o.y = y*s;
  o.z = z*s;
  return o;
};


point point::operator-(point p) {
  point o;
  o.x = this->x - p.x;
  o.y = this->y - p.y;
  o.z = this->z - p.z;
  return o;
};

point point::operator+(point p) {
  point o;
  o.x = x + p.x;
  o.y = y + p.y;
  o.z = z + p.z;
  return o;
};

//--- PARTICLE FUNCTIONS ----
void particle::operator=(point p) {
  x=p.x;
  y=p.y;
  z=p.z;
};

ostream &operator<<(ostream &out, point &p) {
  out << "["<<p.x<<","<<p.y<<","<<p.z<<"]";
  return out;
};

particle::particle() {
  charge=mw=radius=0;
  id=-1;
  hydr=false;
};

spherical::spherical(double a, double b, double c) {
  r=a;
  theta=b;
  phi=c;
};
