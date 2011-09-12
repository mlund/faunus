#include "faunus/point.h"
#include "faunus/slump.h"
#include <faunus/container.h>
#include <faunus/species.h>

namespace Faunus {

  /********************************
    C A R T E S I A N  P O I N T
   ********************************/

  point::point() { clear(); }

  void point::clear() { x=y=z=0; }

  point::point(double xx, double yy, double zz) {
    x=xx;
    y=yy;
    z=zz;
  }

  double point::dot(const point &p) const { return (x*p.x + y*p.y + z*p.z); }

  double point::len() const {
    double l2=x*x+y*y+z*z;
    return (l2!=0) ? sqrt(l2) : 0;
  }

  void point::ranunit(random &ran) {
    point u;
    double r=2;
    while (r > 1.) { //Generate a random unit vector
      u.x=2*ran.random_half();
      u.y=2*ran.random_half();
      u.z=2*ran.random_half();
      r=sqrt(u.x*u.x+u.y*u.y+u.z*u.z);
    }
    x=u.x/r;
    y=u.y/r;
    z=u.z/r;
  }

  point point::operator-() const {
    point o=*this;
    o.x=-x; o.y=-y; o.z=-z;
    return o;
  }

  point & point::operator+=(const point &p) {
    x+=p.x;
    y+=p.y;
    z+=p.z;
    return *this;
  }

  point point::operator+(const point &p) const {
    point o=*this;
    o+=p;
    return o;
  }

  point point::operator-(const point &p) const {
    point o=*this;
    o+=-p;
    //o.x=-o.x;
    //o.y=-o.y;
    //o.z=-o.z;
    return o;
  }

  point& point::operator*=(double s) {
    x*=s;
    y*=s;
    z*=s;
    return *this;
  }

  point point::operator*(double s) const {
    point o=*this;
    o*=s;
    return *this;
  }

  bool point::operator==(const point& p) const {
    return (*this == p);
  }

  string point::str() {
    std::stringstream s;
    s.setf(std::ios::fixed);
    s.precision(2);
    s << x << "," << y << "," << z;
    return "[" + s.str() + "]";
  }

  std::ostream &operator<<(std::ostream &o, const point &p) {
    o << p.x << " " << p.y << " " << p.z;
    return o;
  }

  point & point::operator<<(std::istream &in) {
    in >> x >> y >> z;
    return *this;
  }

  /********************
    P A R T I C L E
   ********************/

  pointparticle::pointparticle() {
    pointparticle::clear();
  }

  void pointparticle::clear() {
    point::clear();
    charge=mw=radius=0;
    hydrophobic=false;
    id=0;
  }

  pointparticle& pointparticle::operator=(const point &p ) {
    x=p.x;
    y=p.y;
    z=p.z;
    return *this;  // return a reference to myself
  }

  double pointparticle::volume() const {
    return (4./3.)*M_PI*radius*radius*radius;
  }

  double pointparticle::mw2vol(double rho) const {
    return 1.6606*rho*mw;
  }

  double pointparticle::mw2rad(double rho) const {
    return pow( mw2vol(rho)*3./4./M_PI, (1/3.) );
  }

  /*!
   * This class tries to deactivate a pointparticle so that certain energy
   * loops does not have to be split. The deactivation is done by
   * \li Setting the charge to zero (no electrostatics)
   * \li Moving the pointparticle *very* far away (p.x=1e9)
   */
  void pointparticle::deactivate() {
    charge=0;
    x=1e9;
  }

  std::ostream &operator<<(std::ostream &o, const pointparticle &p) {
    point b=p;
    o << b << " " << p.charge << " " << p.radius << " " << p.mw << " " << (short)p.id << " " << p.hydrophobic;
    return o;
  }

  pointparticle & pointparticle::operator<<(std::istream &in) {
    short tmp;
    point::operator<<(in);
    in >> charge >> radius >> mw >> tmp >> hydrophobic;
    id = (unsigned char)tmp;
    return *this;
  }

  bool pointparticle::overlap(const pointparticle&a, double r2) const {
    double s=radius+a.radius;
    return (r2<s*s) ? true:false;
  }

  pointparticle pointparticle::operator=(const specdata &d) const {
    pointparticle p;
    p.charge=d.charge;
    p.radius=d.radius;
    p.mw=d.mw;
    p.hydrophobic=d.hydrophobic;
    return p;
  }

  /************
    C I G A R
   ************/

  void cigarparticle::rotate(const geometrybase &c, const point &v, double angle) { //!< Rotate around a vector
  }

  void cigarparticle::translate(const geometrybase &c, const point &v) {             //!< Translate along a vector
  }

  void cigarparticle::scale(const geometrybase &c, double v) {                       //!< Volume scaling
  }

  bool cigarparticle::overlap(const cigarparticle&a, double r2) const {
    return false;
  }

  cigarparticle & cigarparticle::operator<<(std::istream &in) {
    pointparticle::operator<<(in);
    omega.operator<<(in);
    patch.operator<<(in);
    in >> patchangle >> length;
    return *this;
  }

  cigarparticle cigarparticle::operator+(const point &p) const {
    return *this;
  }

  cigarparticle cigarparticle::operator=(const specdata &d) const {
    pointparticle::operator=(d);
    return *this;
  }

  cigarparticle& cigarparticle::operator=(const pointparticle &p ) {
    pointparticle::operator=(p);
    return *this;
  }

  std::ostream &operator<<(std::ostream &o, const cigarparticle &p) {
    o << pointparticle(p)
      << " " << p.omega << " " << p.patch
      << " " << p.patchangle << " " << p.length;
    return o;
  }

}//namespace
