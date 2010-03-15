#include "faunus/point.h"

namespace Faunus {

  /***********************
     H Y P E R P O I N T
   ***********************/

  /*!
   * \param du test
   * \param dv test
   * \param dw test
   */
  void hyperpoint::move(double du, double dv, double dw) {
    double nz1, nz2, nz3, nz4,
           tz1, tz2, tz3, tz4,
           rho=du, omega=dv, fi=dw;
    nz1=sqrt(1.-rho*rho);
    nz2=nz1*cos(fi);
    nz1=nz1*sin(fi);
    nz3=rho*sin(omega);
    nz4=rho*cos(omega);

    hyperpoint e1,e2,e3,te1,te2,te3;
    double fact1,fact2,fact3,nabla_nb,fi_nb;
    nabla_nb=slp.random_one()*2.*acos(-1.);
    fi_nb=acos(slp.random_one());
    e1.z1=cos(nabla_nb);
    e1.z2=sin(nabla_nb);
    e1.z3=0.;
    e1.z4=0.;
    e2.z1=-cos(fi_nb)*sin(nabla_nb);
    e2.z2=cos(fi_nb)*cos(nabla_nb);
    e2.z3=sin(fi_nb);
    e2.z4=0.;
    e3.z1=sin(fi_nb)*sin(nabla_nb);
    e3.z2=-sin(fi_nb)*cos(nabla_nb);
    e3.z3=cos(fi_nb);
    e3.z4=0.;

    // First create a random orthonormal basis set at North Pole
    fact1=e1.z1*z1
      +e1.z2*z2
      +e1.z3*z3;
    te1.z1=e1.z1-1./(1.+z4)*fact1*z1;
    te1.z2=e1.z2-1./(1.+z4)*fact1*z2;
    te1.z3=e1.z3-1./(1.+z4)*fact1*z3;
    te1.z4=e1.z4-1./(1.+z4)*fact1*(z4+1.);
    fact2=e2.z1*z1
      +e2.z2*z2
      +e2.z3*z3;
    te2.z1=e2.z1-1./(1.+z4)*fact2*z1;
    te2.z2=e2.z2-1./(1.+z4)*fact2*z2;
    te2.z3=e2.z3-1./(1.+z4)*fact2*z3;
    te2.z4=e2.z4-1./(1.+z4)*fact2*(z4+1.);
    fact3=e3.z1*z1
      +e3.z2*z2
      +e3.z3*z3;
    te3.z1=e3.z1-1./(1.+z4)*fact3*z1;
    te3.z2=e3.z2-1./(1.+z4)*fact3*z2;
    te3.z3=e3.z3-1./(1.+z4)*fact3*z3;
    te3.z4=e3.z4-1./(1.+z4)*fact3*(z4+1.);

    // Then move it to point of z1,z2,z3,z4
    tz1=nz1*te1.z1+nz2*te2.z1+nz3*te3.z1+nz4*z1;
    tz2=nz1*te1.z2+nz2*te2.z2+nz3*te3.z2+nz4*z2;
    tz3=nz1*te1.z3+nz2*te2.z3+nz3*te3.z3+nz4*z3;
    tz4=nz1*te1.z4+nz2*te2.z4+nz3*te3.z4+nz4*z4;

    // Update the point
    z1=tz1;
    z2=tz2;
    z3=tz3;
    z4=tz4;
  }

  std::ostream & operator<<(std::ostream &o, hyperpoint &p) {
    o << p.z1 << " " << p.z2 << " " << p.z3 << " " << p.z4;
    return o;
  }

  hyperpoint & hyperpoint::operator<<(std::istream &in) {
    in >> z1 >> z2 >> z3 >> z4;
    return *this;
  }

  /********************************
    C A R T E S I A N  P O I N T
   ********************************/

  point::point() {
    clear();
  }

  point::point(double xx, double yy, double zz) {
    x=xx;
    y=yy;
    z=zz;
  }

  void point::clear() {
    x=y=z=0;
  }

  double point::dot(const point &p) const {
    return (x*p.x + y*p.y + z*p.z);
  }

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

  point & point::operator+=(const point &p) {
    x += p.x;
    y += p.y;
    z += p.z;
    return *this;
  }

  point point::operator-() {
    point o;
    o.x=-x; o.y=-y; o.z=-z;
    return o;
  }

  point point::operator*(const point p) {
    point o;
    o.x = p.x * x;
    o.y = p.y * y;
    o.z = p.z * z;
    return o;
  }

  point point::operator*(double s) const {
    point o;
    o.x = x*s;
    o.y = y*s;
    o.z = z*s;
    return o;
  }

  point point::operator+(double d) {
    point o;
    o.x = x+d;
    o.y = y+d;
    o.z = z+d;
    return o;
  }

  point point::operator-(const point p) const {
    point o;
    o.x = this->x - p.x;
    o.y = this->y - p.y;
    o.z = this->z - p.z;
    return o;
  }

  point point::operator+(const point p) {
    point o;
    o.x = x + p.x;
    o.y = y + p.y;
    o.z = z + p.z;
    return o;
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

  std::ostream &operator<<(std::ostream &o, point &p) {
#ifdef HYPERSPHERE
    hyperpoint hp=p;
    o << hp << " ";
#endif
    o << p.x << " " << p.y << " " << p.z;
    return o;
  }

  point & point::operator<<(std::istream &in) {
#ifdef HYPERSPHERE
    hyperpoint::operator<<(in);
#endif
    in >> x >> y >> z;
    return *this;
  }

  /********************
    P A R T I C L E
   ********************/

  particle::particle() {
    particle::clear();
  }

  void particle::clear() {
    point::clear();
    charge=mw=radius=0;
    hydrophobic=false;
    id=0;
  }

  particle& particle::operator=(const point& p) {
    x=p.x;
    y=p.y;
    z=p.z;
    return *this;  // return a reference to myself
  }

  double particle::volume() const {
    return (4./3.)*M_PI*radius*radius*radius;
  }

  double particle::mw2vol(double rho) const {
    return 1.6606*rho*mw;
  }

  double particle::mw2rad(double rho) const {
    return pow( mw2vol(rho)*3./4./M_PI, (1/3.) );
  }

  /*!
   * This class tries to deactivate a particle so that certain energy
   * loops does not have to be split. The deactivation is done by
   * \li Setting the charge to zero (no electrostatics)
   * \li Moving the particle *very* far away (p.x=1e9)
   * \li Moving hyperpoints to the opposide side of the sphere (no overlap, minimum vdW)
   */
  void particle::deactivate() {
    charge=0;
#ifndef HYPERSPHERE
    x=1e9;
#else
    z1*=-0.9999;
    z2*=-0.9999;
    z3*=-0.9999;
    z4*=-0.9999;
#endif
  }

  std::ostream &operator<<(std::ostream &o, particle &p) {
    point b=p;
    o << b << " " << p.charge << " " << p.radius << " " << p.mw << " " << p.id << " " << p.hydrophobic;
    return o;
  }

  particle & particle::operator<<(std::istream &in) {
    point::operator<<(in);
    in >> charge >> radius >> mw >> id >> hydrophobic;
    return *this;
  }

  //---------------- SPHERICAL -----------------
  spherical::spherical(double radial, double zenith, double azimuthal) {
    r=radial;
    theta=zenith;
    phi=azimuthal;
  }
}//namespace
