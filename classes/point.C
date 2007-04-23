#include "point.h"

//Constructor
point::point() { clear(); }

void point::clear() { x=y=z=0; };
double point::dot(point &p) { return (x*p.x + y*p.y + z*p.z); }

//! Length of vector
double point::len() {
  double l2=x*x+y*y+z*z;
  return (l2!=0) ? sqrt(l2) : 0;
}

void point::operator+=(point p) {
  x += p.x;
  y += p.y;
  z += p.z;
}

// Sign reversal
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

point point::operator+(double d) {
  point o;
  o.x = x+d;
  o.y = y+d;
  o.z = z+d;
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

/*
void point::random_unitv() {
  double r=2;
  while (r > 1.) { //random unit vector
    x=random_one();
    y=random_one();
    z=random_one();
    r=sqrt(x*x+y*y+z*z);
  };
  x=x/r;
  y=y/r;
  z=z/r;
};
*/

//--- PARTICLE FUNCTIONS ----
void particle::operator=(point p) {
  x=p.x;
  y=p.y;
  z=p.z;
};

string point::str() {
  stringstream s;
  s.setf(ios::fixed);
  s.precision(2);
  s << x << "," << y << "," << z;
  return "[" + s.str() + "]";
};

ostream &operator<<(ostream &out, point &p) {
  out << p.str();
  return out;
};

particle::particle() {
  charge=mw=radius=0;
  id=GHOST;
};

spherical::spherical(double radial, double zenith, double azimuthal) {
  r=radial;
  theta=zenith;
  phi=azimuthal;
};
