#include <faunus/geometry.h>
#include <faunus/inputfile.h>
#include <faunus/species.h>
#include <faunus/physconst.h>
#include <faunus/textio.h>
#include <faunus/group.h>
#include <Eigen/Core>
#include <Eigen/Geometry>

namespace Faunus {
  namespace Geometry {

    using namespace Faunus::textio;

    Geometrybase::~Geometrybase() {}

    double Geometrybase::dist(const Point &p1, const Point &p2) const {
      return sqrt(sqdist(p1,p2));
    }

    void Geometrybase::scale(Point &a, Point &s, const double xyz, const double xy) const {
      assert(!"Scaling function unimplemented for this geometry");
    }

    string Geometrybase::info(char w) {
      using namespace textio;
      std::ostringstream o;
      o << pad(SUB,w, "Boundary") << name << endl
        << pad(SUB,w, "Volume") << getVolume() << _angstrom << cubed
        << " = " << getVolume()/1e3 << " nm" << cubed
        << " = " << getVolume()/1e27 << " liters" << endl
        << _info(w);
      return o.str();
    }

    void Geometrybase::setVolume(double volume) {
      assert( volume>0 && "Zero volume not allowed!");
      _setVolume( volume );
      assert( std::abs( (volume-getVolume())/volume )<1e-9
          && "setVolume() and/or getVolume() is broken!" );
    }

    double Geometrybase::getVolume() const {
      return _getVolume();
    }

    bool Geometrybase::save(string file) {
      std::ofstream fout( file.c_str());
      if (fout) {
        fout.precision(10);
        fout <<  getVolume() << endl;
        return true;
      }
      return false;
    }

    /*!
     * \param file Filename
     * \param resize True if the current geometry should be resized to match file content (default: false)
     */
    bool Geometrybase::load(string file, bool resize) {
      std::ifstream f( file.c_str() );
      if (f) {
        double vol;
        f >> vol;
        setVolume(vol);
        return true;
      }
      std::cerr << "!! Geometry data NOT read from file " << file << endl;
      return false;
    }

    Sphere::Sphere(double radius) {
      setRadius(radius);
      len = Point(r,diameter,0);
    }

    /**
     * The InputMap is scanned for the following parameters:
     *
     * Key               | Description
     * :---------------- | :-------------------------------------------------------
     * `sphere_radius`   | Sphere radius [angstrom]
     */
    Sphere::Sphere(InputMap &in, string prefix)  {
      setRadius(
          in.get<double>(prefix+"_radius", -1.0, "Spherical container radius (A)") );
    }

    void Sphere::setRadius(double radius) {
      assert(radius>0 && "Radius must be larger than zero.");
      name="Spherical";
      r = radius; 
      r2 = r*r; 
      diameter = 2*r; 
      len = Point(r,diameter,0);
    }

    double Sphere::_getVolume() const {
      return 4/3.*pc::pi*std::pow(r,3);
    }

    void Sphere::_setVolume(double vol) {
      setRadius( cbrt( 3*vol/(4*pc::pi) ) );
    }

    bool Sphere::setlen(const Point &l) {
      if (l.x()<=0) return false;
      Sphere::setRadius(l.x());
      return true;
    }

    void Sphere::scale(Point &a, Point &s, const double xyz=1, const double xy=1) const {
      assert( getVolume()>0 );
      double newradius = cbrt( 3*xyz*xyz*xyz*getVolume()/(4*pc::pi) );
      a = a * (newradius/r);
    }

    string Sphere::_info(char w) {
      std::ostringstream o;
      o << pad(SUB,w,"Radius") << r << textio::_angstrom << endl;
      return o.str();
    }

    void Sphere::randompos(Point &a) {
      do {
        a.x() = (slp()-0.5)*diameter;
        a.y() = (slp()-0.5)*diameter;
        a.z() = (slp()-0.5)*diameter;
      } while ( a.squaredNorm()>r2 );
    }

    bool Sphere::collision(const Point &a, double radius, collisiontype type) const {
      return (a.squaredNorm()>r2) ? true:false;
    }

    /**
     * The InputMap is scanned for the following parameters:
     *
     * Key               | Description
     * :---------------- | :-------------------------------------------------------
     * `cuboid_len`      | Uniform sidelength [angstrom]. If negative, continue to...
     * `cuboid_xlen`     | x sidelength [angstrom]
     * `cuboid_ylen`     | y sidelength [angstrom]
     * `cuboid_zlen`     | z sidelength [angstrom]
     * `cuboid_scaledir` | Isobaric scaling directions (`XYZ`=isotropic, `XY`=xy only).
     */
    Cuboid::Cuboid(InputMap &in) {
      name="Cuboid";
      string scaledirstr = in.get<string>("cuboid_scaledir","XYZ");
      if (scaledirstr=="XY")
        scaledir=XY;
      else
        scaledir=XYZ;
      double cubelen=in.get<double>("cuboid_len",-1, name+" sidelength (AA)");
      if (cubelen<1e-6) {
        len.x()=in.get<double>("cuboid_xlen",0);
        len.y()=in.get<double>("cuboid_ylen",0);
        len.z()=in.get<double>("cuboid_zlen",0);
      } else
        len.x()=len.y()=len.z()=cubelen;
      setlen(len);
    }

    bool Cuboid::setlen(const Point &l) {
      //if ( l.squaredNorm() < 1e-6 ) {
      //  std::cerr << "Error: cuboid volume is zero. " << endl;
      //  exit(1);
      //}
      assert(l.x()>0 && l.y()>0 && l.z()>0);
      if (l.x()<=0||l.y()<=0||l.z()<=0) 
        return false;
      len = l;                    // Cuboid sidelength
      len_half=l*0.5;             // half Cuboid sidelength
      len_inv.x()=1/len.x();      // inverse Cuboid side length
      len_inv.y()=1/len.y();
      len_inv.z()=1/len.z();
      return true;
    }

    void Cuboid::_setVolume(double newV) {
      double xyz=1.0, xy=1.0;
      Point s = Point(1,1,1);
      if (scaledir==XYZ) 
        xyz = cbrt(newV/getVolume());
      if (scaledir==XY) 
        xy = sqrt(newV/getVolume());
      scale(len,s,xyz,xy);
      setlen(len);
    }

    double Cuboid::_getVolume() const {
      return len.x()*len.y()*len.z();
    }

    string Cuboid::_info(char w) {
      std::ostringstream o;
      o << pad(SUB,w, "Sidelengths")
        << len.transpose() << " ("+textio::angstrom+")" << endl
        << pad(SUB,w, "Scale directions") << scaledirstr << endl;
      return o.str();
    }

    Point Cuboid::randompos() {
      Point m;
      randompos(m);
      return m;
    }

    void Cuboid::randompos(Point &m) {
      m.x() = slp.half()*len.x();
      m.y() = slp.half()*len.y();
      m.z() = slp.half()*len.z();
    }

    bool Cuboid::save(string file) {
      if ( Geometrybase::save(file) ) {
        std::ofstream fout( file.c_str(), std::ios_base::app);
        if (fout) {
          fout.precision(10);
          fout << len.transpose() << endl;
          return true;
        }
      }
      return false;
    }

    bool Cuboid::load(string file, bool resize) {
      Point l;
      if ( Geometrybase::load(file, resize) ) {
        std::ifstream f( file.c_str() );
        if (f) {
          f >> l.x() >> l.y() >> l.z();
          setlen(l);
          return true;
        }
      }
      std::cerr << "# Container data NOT read from file " << file << endl;
      return false;
    }

    void Cuboid::scale(Point &a, Point &s, const double xyz=1, const double xy=1) const {
      assert( getVolume()>0 );
      if (scaledir==XYZ)
        a = Point(a.x()*s.x()*xyz,a.y()*s.y()*xyz,a.z()*s.z()*xyz);
      if (scaledir==XY) 
        a = Point(a.x()*s.x()*xy,a.y()*s.y()*xy,a.z()*s.z());
    }

    Cuboidslit::Cuboidslit(InputMap &in) : Cuboid(in) {
      name="Cuboid XY-periodicity";
    }

    /**
     * @param length Length of the Cylinder (angstrom)
     * @param radius Radius of the Cylinder (angstrom)
     */
    Cylinder::Cylinder(double length, double radius) {
      init(length, radius);
    }

    /**
     * The InputMap is scanned for the following parameters:
     *
     * Key               | Description
     * :---------------- | :-------------------------------------------------------
     * `cylinder_len`    | Cylinder length [angstrom]
     * `cylinder_radius` | Cylinder radius [angstrom] 
     */
    Cylinder::Cylinder(InputMap &in) {
      init(
          in.get<double>("cylinder_len", 0), in.get<double>("cylinder_radius", 0) );
    }

    void Cylinder::init(double length, double radius) {
      name="Cylindrical (hard ends)";
      assert(length>0 && radius>0 && "Cylinder length and radius must be bigger than zero.");
      _len=length;
      setVolume( pc::pi*radius*radius*_len );
    }

    /**
     * Dummy function not be used other than for compatibility
     * reasons. Sets length to x component of vector.
     */
    bool Cylinder::setlen(const Point &l) {
      init( l.x(), r); 
      return true;
    }

    void Cylinder::_setVolume(double newV) {
      r2=newV/(pc::pi*_len);
      r=sqrt(r2);
      diameter=2*r;
      _halflen=_len/2;
    }

    double Cylinder::_getVolume() const {
      return r2*pc::pi*_len;
    }

    void Cylinder::boundary(Point &p) const {}

    void Cylinder::randompos(Point &m) {
      double l=r2+1;
      m.z() = slp.half()*_len;
      while (l>r2) {
        m.x() = slp.half()*diameter;
        m.y() = slp.half()*diameter;
        l=m.x()*m.x()+m.y()*m.y();
      }
    }

    bool Cylinder::collision(const Point &a, double radius, collisiontype type) const {
      assert( (_halflen-_len/2)<1e-9 && "Cylinder length initialization problems" );
      if ( a.z()<-_halflen ) return true;
      if ( a.z()>_halflen ) return true;
      if ( a.x()*a.x()+a.y()*a.y()>r2 ) return true;
      return false;
    }

    string Cylinder::_info(char w) {
      std::ostringstream o;
      o << pad(SUB,w, "Length") << _halflen*2 << textio::_angstrom << endl
        << pad(SUB,w, "Radius") << sqrt(r2) << textio::_angstrom << endl;
      return o.str();
    }

    PeriodicCylinder::PeriodicCylinder(double length, double radius) : Cylinder(length,radius) {
      name="Cylindrical (periodic ends)";
    }

    PeriodicCylinder::PeriodicCylinder(InputMap &in) : Cylinder(in) {
      name="Cylindrical (periodic ends)";
    }

    void PeriodicCylinder::boundary(Point &a) const {
      if (std::abs(a.z())>_halflen)
        a.z()-=_len*anint(a.z()/_len);
    }

  }//namespace geometry
}//namespace
