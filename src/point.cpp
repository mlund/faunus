#include "faunus/point.h"
#include "faunus/slump.h"
#include <faunus/geometry.h>
#include <faunus/species.h>
#include <faunus/physconst.h>

namespace Faunus {

  /********************************
    C A R T E S I A N  P O I N T
   ********************************/

  /*!
   * Upon construction x,y,z are set to zero
   */
  Point::Point() : x(0), y(0), z(0) {}

  Point::Point(double xx, double yy, double zz) : x(xx), y(yy), z(zz) {}

  void Point::clear() { x=y=z=0; }

  double Point::dot(const Point &p) const { return (x*p.x + y*p.y + z*p.z); }

  double Point::len() const {
    double l2=x*x+y*y+z*z;
    return (l2!=0) ? sqrt(l2) : 0;
  }

  void Point::ranunit(RandomBase &ran) {
    Point u;
    double r=2;
    while (r > 1.) { //Generate a random unit vector
      u.x=2*ran.randHalf();
      u.y=2*ran.randHalf();
      u.z=2*ran.randHalf();
      r=sqrt(u.x*u.x+u.y*u.y+u.z*u.z);
    }
    x=u.x/r;
    y=u.y/r;
    z=u.z/r;
  }

  Point Point::operator-() const {
    Point o=*this;
    o.x=-x; o.y=-y; o.z=-z;
    return o;
  }

  Point & Point::operator+=(const Point &p) {
    x+=p.x;
    y+=p.y;
    z+=p.z;
    return *this;
  }

  Point Point::operator+(const Point &p) const {
    Point o=*this;
    o+=p;
    return o;
  }

  Point Point::operator-(const Point &p) const {
    Point o=*this;
    o+=-p;
    //o.x=-o.x;
    //o.y=-o.y;
    //o.z=-o.z;
    return o;
  }

  Point& Point::operator*=(double s) {
    x*=s;
    y*=s;
    z*=s;
    return *this;
  }

  const Point Point::operator*(double s) const {
    Point o=*this;
    o*=s;
    return o;
  }

  bool Point::operator==(const Point& p) const {
    if (&p==&(*this))
      return true;
    if (p.x!=x) return false;
    if (p.y!=y) return false;
    if (p.z!=z) return false;
    return true;
  }

  string Point::str() {
    std::stringstream s;
    s.setf(std::ios::fixed);
    s.precision(2);
    s << x << "," << y << "," << z;
    return "[" + s.str() + "]";
  }

  std::ostream &operator<<(std::ostream &o, const Point &p) {
    o << p.x << " " << p.y << " " << p.z;
    return o;
  }

  Point & Point::operator<<(std::istream &in) {
    in >> x >> y >> z;
    return *this;
  }

  /********************
    P A R T I C L E
   ********************/

  /*!
   * Upon construction data is zeroed.
   */
  PointParticle::PointParticle() : charge(0), radius(0), mw(0), id(0), hydrophobic(false) {}

  void PointParticle::clear() {
    Point::clear();
    charge=mw=radius=0;
    hydrophobic=false;
    id=0;
  }

  PointParticle& PointParticle::operator=(const Point &p) {
    x=p.x;
    y=p.y;
    z=p.z;
    return *this;  // return a reference to myself
  }

  double PointParticle::volume() const {
    return (4./3.)*pc::pi*radius*radius*radius;
  }

  double PointParticle::mw2vol(double rho) const {
    return 1.6606*rho*mw;
  }

  double PointParticle::mw2rad(double rho) const {
    return pow( mw2vol(rho)*3./4./pc::pi, (1/3.) );
  }

  /*!
   * This class tries to deactivate a pointparticle so that certain energy
   * loops does not have to be split. The deactivation is done by
   * \li Setting the charge to zero (no electrostatics)
   * \li Moving the pointparticle *very* far away (p.x=1e9)
   */
  void PointParticle::deactivate() {
    charge=0;
    x=1e9;
  }

  std::ostream &operator<<(std::ostream &o, const PointParticle &p) {
    Point b=p;
    o << b << " " << p.charge << " " << p.radius << " " << p.mw << " " << (PointParticle::Tid)p.id << " " << p.hydrophobic;
    return o;
  }

  PointParticle & PointParticle::operator<<(std::istream &in) {
    int tmp; // avoid char text output
    Point::operator<<(in);
    in >> charge >> radius >> mw >> tmp >> hydrophobic;
    id = (Tid)tmp;
    return *this;
  }

  bool PointParticle::overlap(const PointParticle&a, double r2) const {
    double s=radius+a.radius;
    return (r2<s*s) ? true:false;
  }

  PointParticle& PointParticle::operator=(const AtomData &d) {
    id=d.id;
    charge=d.charge;
    radius=d.radius;
    mw=d.mw;
    hydrophobic=d.hydrophobic;
    return *this;
  }

  /*****************************
    S P H E R O C Y L I N D E R
   *****************************/

  void CigarParticle::rotate(const Geometrybase &c, const Point &v, double angle) { //!< Rotate around a vector
  }

  void CigarParticle::translate(const Geometrybase &c, const Point &v) {             //!< Translate along a vector
  }

  void CigarParticle::scale(const Geometrybase &c, double v) {                       //!< Volume scaling
  }

  bool CigarParticle::overlap(const CigarParticle&a, double r2) const {
    return false;
  }

  CigarParticle & CigarParticle::operator<<(std::istream &in) {
    PointParticle::operator<<(in);
    omega.operator<<(in);
    patch.operator<<(in);
    in >> patchangle >> length;
    return *this;
  }

  CigarParticle CigarParticle::operator+(const Point &p) const {
    return *this;
  }

  CigarParticle& CigarParticle::operator=(const Point &p) {
    PointParticle::operator=(p);
    return *this;
  }

  CigarParticle& CigarParticle::operator=(const AtomData &d) {
    PointParticle::operator=(d);
    return *this;
  }

  CigarParticle& CigarParticle::operator=(const PointParticle &p) {
    PointParticle::operator=(p);
    return *this;
  }

  std::ostream &operator<<(std::ostream &o, const CigarParticle &p) {
    o << PointParticle(p)
      << " " << p.omega << " " << p.patch
      << " " << p.patchangle << " " << p.length;
    return o;
  }

}//namespace
