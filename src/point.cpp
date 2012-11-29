#include <faunus/point.h>
#include <faunus/slump.h>
#include <faunus/geometry.h>
#include <faunus/species.h>
#include <faunus/physconst.h>

namespace Faunus {

  /********************************
    C A R T E S I A N  P O I N T
   ********************************/
   
  void Point::clear() { setZero(); }

  Point::Tcoord Point::len() const {
    auto r2 = dot(*this);
    assert( r2==squaredNorm() );
    return (r2>0) ? std::sqrt(r2) : 0;
  }

  /*!
   * This uses the von Neumann method described in Allen and Tildesley page 349.
   */
  void Point::ranunit(RandomBase &ran) {
    Point u;
    Tcoord r2;
    do {
      u.x()=2*ran.randHalf();
      u.y()=2*ran.randHalf();
      u.z()=2*ran.randHalf();
      r2=u.squaredNorm();
    } while (r2>1);
    *this = u/std::sqrt(r2);
    assert(std::abs(this->len()-1)<1e-7); // is it really a unit vector?
  }

  Point & Point::operator<<(std::istream &in) {
    in >> x() >> y() >> z();
    return *this;
  }

  /*!
   * \param vrot Rotation class where the axis, angle and geometry are expected to be set
   */
  void Point::rotate(Geometry::VectorRotate &vrot) {
    *this = vrot.rotate(*this);
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
    x()=y()=z()=pc::infty;
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
    rot.rotate(dir);
    rot.rotate(patchdir);
    rot.rotate(patchsides[0]);
    rot.rotate(patchsides[1]);
    rot.rotate(chdir);
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
      << " " << p.dir << " " << p.patchdir;
    return o;
  }

  // remember to read ALL data members from stream - and in same order as written!
  CigarParticle& CigarParticle::operator<<(std::istream &in) {
    PointParticle::operator<<(in);
    dir.operator<<(in);
    patchdir.operator<<(in);
    return *this;
  }

}//namespace
