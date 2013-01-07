#include <faunus/point.h>
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

  Point & Point::operator<<(std::istream &in) {
    in >> x() >> y() >> z();
    return *this;
  }

  /*!
   * \param rotator Functor that rotates a point and returns the rotated Point
   *
   * The functor should take care of simulation boundaries (if any) and typically one
   * would want to pass the Geometry::VectorRotate class as in the following example:
   * \code
   * Point a(1,0,0);
   * VectorRotate rotator;
   * rotator.setAxis(geometry, Point(0,0,0), Point(0,0,1), 3.14 ); // rotate pi around 0,0,1
   * a.rotate(rotator);
   * \endcode
   */
  void Point::rotate(RotFunctor rotator) {
    *this = rotator(*this);
  }
  
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

  /**
   * This class tries to deactivate a pointparticle so that certain energy
   * loops does not have to be split. The deactivation is done by
   *
   * - Setting the charge to zero (no electrostatics)
   * - Setting the hydrophobic tag to false
   * - Moving the pointparticle *very* far away
   */
  void PointParticle::deactivate() {
    hydrophobic=false;
    charge=0;
    x()=y()=z()=pc::infty;
  }

  /**
   * This will write all data to given stream. Note that the particle id is converted
   * to a short integer since char output (Tid=char) does not print well on screen. Derived
   * classes should expand on this so that *all* data is written.
   */
  std::ostream& operator<<(std::ostream &o, const PointParticle &p) {
    o << Point(p)
      << " " << p.charge << " " << p.radius << " " << p.mw << " " << (short)p.id << " " << p.hydrophobic;
    return o;
  }

  /**
   * This will read all data from stream in the same order as writte.
   * Note that a short integer is expected for the particle id
   * since chars (Tid=char) does not print well on screen.
   * Derived classes should expand on this so
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

}//namespace
