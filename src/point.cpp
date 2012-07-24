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
   * Data is NOT zeroed upon construction!
   */
  Point::Point() {}

  Point::Point(Tcoord xx, Tcoord yy, Tcoord zz) : x(xx), y(yy), z(zz) {}

  Point::~Point() {}

  void Point::clear() { x=y=z=0; }

  Point::Tcoord Point::dot(const Point &p) const { return (x*p.x + y*p.y + z*p.z); }

  Point::Tcoord Point::len() const {
    auto r2 = dot(*this);
    return (r2>0) ? sqrt(r2) : 0;
  }

  void Point::ranunit(RandomBase &ran) {
    Point u;
    Tcoord r=2;
    while (r>1) { //Generate a random unit vector
      u.x=2*ran.randHalf();
      u.y=2*ran.randHalf();
      u.z=2*ran.randHalf();
      r=sqrt(u.x*u.x+u.y*u.y+u.z*u.z);
    }
    *this = u*(1/r);
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
    return o;
  }

  Point& Point::operator*=(Tcoord s) {
    x*=s;
    y*=s;
    z*=s;
    return *this;
  }

  const Point Point::operator*(Tcoord s) const {
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

  std::ostream &operator<<(std::ostream &o, const Point &p) {
    o << p.x << " " << p.y << " " << p.z;
    return o;
  }

  Point & Point::operator<<(std::istream &in) {
    in >> x >> y >> z;
    return *this;
  }

  /*!
   * \param rot Rotation class where the axis, angle and geometry are expected to be set
   */
  void Point::rotate(Geometry::VectorRotate &vrot) {
    vrot.rotate(*this);
  }

  /*!
   * \param geo Geomtry to use so that boundary conditions can be respected
   * \param a Vector to translate with
   */
  void Point::translate(const Geometry::Geometrybase &geo, const Point &a) {
    assert(&geo!=nullptr);
    (*this)+=a;
    geo.boundary(*this);
  }

  /*!
   * This will perform a volume scaling of the Point by following the algorithm
   * specified in the Geometrybase. Derived classes for more complex particles
   * may override this function.
   */
  void Point::scale(const Geometry::Geometrybase &geo, double newvol) {
    geo.scale(*this, newvol);
  }

  /*
  void Point::vectoriz(vector<double> &v, vectorizeFmt fmt) const {
    v.push_back(x);
    v.push_back(y);
    v.push_back(z);
  }*/

  /********************
    P A R T I C L E
   ********************/

  /*!
   * Upon construction, data is zeroed.
   */
  PointParticle::PointParticle() {
    clear();
  }

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

  /*!
   * This class tries to deactivate a pointparticle so that certain energy
   * loops does not have to be split. The deactivation is done by
   * \li Setting the charge to zero (no electrostatics)
   * \li Setting the hydrophobic tag to false
   * \li Moving the pointparticle *very* far away
   */
  void PointParticle::deactivate() {
    hydrophobic=false;
    charge=0;
    x=y=z=pc::infty;
  }

  /*!
   * This will write all data to given stream. Note that the particle id is converted
   * to a short integer since char output (Tid=char) does not print well on screen. Derived
   * classes should expand on this so that *all* data is written.
   */
  std::ostream& operator<<(std::ostream &o, const PointParticle &p) {
    o << Point(p)
      << " " << p.charge << " " << p.radius << " " << p.mw << " " << (short)p.id << " " << p.hydrophobic;
    return o;
  }

  /*!
   * This will read all data from stream in the same order as writte.
   * Note that a short integer is expected for the particle id
   * since chars (Tid=char) does not print well on screen. Derived classes should expand on this so
   * that *all* data is read.
   */
  PointParticle& PointParticle::operator<<(std::istream &in) {
    short tmp; // avoid char output in readable text files
    Point::operator<<(in);
    in >> charge >> radius >> mw >> tmp >> hydrophobic;
    id = (Tid)tmp;
    return *this;
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

  void CigarParticle::rotate(Geometry::VectorRotate &rot) {
    assert(!"Unimplemented");
  }

  void CigarParticle::translate(const Geometry::Geometrybase &geo, const Point &a) {
    assert(!"Unimplemented");
  }

  void CigarParticle::scale(const Geometry::Geometrybase &geo, double newvol) {
    assert(!"Unimplemented");
  }

  CigarParticle CigarParticle::operator+(const Point &p) const {
    assert(!"Unimplemented");
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

  // remember to write ALL data members to stream!
  std::ostream &operator<<(std::ostream &o, const CigarParticle &p) {
    o << PointParticle(p)
      << " " << p.omega << " " << p.patch
      << " " << p.patchangle << " " << p.length;
    return o;
  }

  // remember to read ALL data members from stream - and in same order as written!
  CigarParticle& CigarParticle::operator<<(std::istream &in) {
    PointParticle::operator<<(in);
    omega.operator<<(in);
    patch.operator<<(in);
    in >> patchangle >> length;
    return *this;
  }

}//namespace
