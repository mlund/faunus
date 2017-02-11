#include <faunus/geometry.h>
#include <faunus/inputfile.h>
#include <faunus/species.h>
#include <faunus/physconst.h>
#include <faunus/textio.h>
#include <faunus/group.h>
#include <Eigen/Core>
#include <Eigen/Geometry>

namespace Faunus
{
  namespace Geometry
  {

    using namespace Faunus::textio;

    Geometrybase::Geometrybase() {}

    Geometrybase::Geometrybase( const string &name ) : name(name)
    {
    }

    Geometrybase::~Geometrybase() {}

    double Geometrybase::dist( const Point &p1, const Point &p2 ) const
    {
        return sqrt(sqdist(p1, p2));
    }

    void Geometrybase::scale( Point &a, Point &s, const double xyz, const double xy ) const
    {
        assert(!"Scaling function unimplemented for this geometry");
    }

    Cuboid Geometrybase::inscribe() const
    {
        assert(!"Inscribe function not implemented for this geometry");
        return Cuboid();
    }

    string Geometrybase::info( char w )
    {
        using namespace textio;
        std::ostringstream o;
        o << pad(SUB, w, "Boundary") << name << endl
          << pad(SUB, w, "Volume") << getVolume() << _angstrom << cubed
          << " = " << getVolume() / 1e3 << " nm" << cubed
          << " = " << getVolume() / 1e27 << " liters" << endl
          << _info(w);
        return o.str();
    }

    void Geometrybase::setVolume( double volume )
    {
        _setVolume(volume);
    }

    double Geometrybase::getVolume(int dim) const
    {
        assert(dim==1 || dim==2 || dim==3);
        return _getVolume(dim);
    }

    Sphere::Sphere( double radius ) : Geometrybase("Sphere")
    {
        len = Point(r, diameter, 0);
    }

    /**
     * Key        | Description
     * :--------- | :-----------------------
     * `radius`   | Sphere radius [angstrom]
     */
    Sphere::Sphere( Tmjson &j ) : Geometrybase("Sphere")
      {
          try {
              setRadius( j.at("radius") );
          }
          catch(std::exception& e) {
              throw std::runtime_error(name + ": " + e.what());
          }
      }


    void Sphere::setRadius( double radius )
    {
        assert(radius > 0 && "Radius must be larger than zero.");
        r = radius;
        r2 = r * r;
        diameter = 2 * r;
        len = Point(r, diameter, 0);
    }

    /**
     * - 3D: \f$ 4\pi r^3/3 \f$ (default)
     * - 2D: \f$ \pi r^2 \f$
     * - 1D: \f$ 2r \f$
     */
    double Sphere::_getVolume(int dim) const
    {
        if (dim==1)
            return 2*r;
        if (dim==2)
            return pc::pi * r * r;
        return 4 / 3. * pc::pi * std::pow(r, 3);
    }

    void Sphere::_setVolume( double vol )
    {
        setRadius(cbrt(3 * vol / (4 * pc::pi)));
    }

    void Sphere::setlen( const Point &l )
    {
        Sphere::setRadius(l.x());
        if ( getVolume() < 1e-9 )
            throw std::runtime_error("Sphere volume is zero.");
    }

    void Sphere::scale( Point &a, Point &s, const double xyz = 1, const double xy = 1 ) const
    {
        assert(getVolume() > 0);
        double newradius = cbrt(3 * xyz * xyz * xyz * getVolume() / (4 * pc::pi));
        a = a * (newradius / r);
    }

    string Sphere::_info( char w )
    {
        std::ostringstream o;
        o << pad(SUB, w, "Radius") << r << textio::_angstrom << endl;
        return o.str();
    }

    void Sphere::randompos( Point &a )
    {
        do
        {
            a.x() = (slump() - 0.5) * diameter;
            a.y() = (slump() - 0.5) * diameter;
            a.z() = (slump() - 0.5) * diameter;
        }
        while ( a.squaredNorm() > r2 );
    }

    bool Sphere::collision( const Point &a, double radius, collisiontype type ) const
    {
        return (a.squaredNorm() > r2) ? true : false;
    }

    Cuboid Sphere::inscribe() const
    {
        Cuboid c;
        c.setlen({diameter, diameter, diameter});
        return c;
    }

    Cuboid Cuboid::inscribe() const
    {
        return *this;
    }

    Cuboid::Cuboid() {}

    void Cuboid::setlen( const Point &l )
    {
        len = l;                    // Cuboid sidelength
        len_half = l * 0.5;             // half Cuboid sidelength
        len_inv.x() = 1 / len.x();      // inverse Cuboid side length
        len_inv.y() = 1 / len.y();
        len_inv.z() = 1 / len.z();
        if ( getVolume() < 1e-9 )
            throw std::runtime_error("Cuboid volume is zero");
    }

    void Cuboid::_setVolume( double newV )
    {
        double xyz = 1.0, xy = 1.0;
        Point s = Point(1, 1, 1);
        if ( scaledir == XYZ )
            xyz = cbrt(newV / getVolume());
        if ( scaledir == XY )
            xy = sqrt(newV / getVolume());
        scale(len, s, xyz, xy);
        setlen(len);
    }

    /**
     * - 3D: \f$ l_x l_y l_z \f$ (default)
     * - 2D: \f$ l_x l_y \f$
     * - 1D: \f$ l_z \f$
     */
    double Cuboid::_getVolume(int dim) const
    {
        if (dim==1)
            return len.z();
        if (dim==2)
            return len.x() * len.y();
        return len.x() * len.y() * len.z();
    }

    string Cuboid::_info( char w )
    {
        std::ostringstream o;
        o << pad(SUB, w, "Sidelengths")
          << len.transpose() << " (" + textio::angstrom + ")" << endl
          << pad(SUB, w, "Scale directions") << scaledirstr << endl;
        return o.str();
    }

    Point Cuboid::randompos()
    {
        Point m;
        randompos(m);
        return m;
    }

    void Cuboid::randompos( Point &m )
    {
        m.x() = slump.half() * len.x();
        m.y() = slump.half() * len.y();
        m.z() = slump.half() * len.z();
    }

    void Cuboid::scale( Point &a, Point &s, const double xyz = 1, const double xy = 1 ) const
    {
        assert(getVolume() > 0);
        if ( scaledir == XYZ )
            a = Point(a.x() * s.x() * xyz, a.y() * s.y() * xyz, a.z() * s.z() * xyz);
        if ( scaledir == XY )
            a = Point(a.x() * s.x() * xy, a.y() * s.y() * xy, a.z() * s.z());
    }

    /**
     * The json object is scanned for the following parameters:
     *
     * Key           | Description
     * :------------ | :-----------------------------------------------------------
     * `length`      | Array of sidelengths _or_ single length for cube [angstrom]
     * `scaledir`    | Isobaric scaling directions (`XYZ`=isotropic, `XY`=xy only).
     */
    Cuboid::Cuboid( Tmjson &j ) : Geometrybase("Cuboid")
    {
        try {
            if (j.is_object()) {
                scaledirstr = j.value<string>("scaledir", "XYZ");
                scaledir = (scaledirstr == "XY") ? XY : XYZ;
                auto m = j.at("length");
                if (m.is_number()) {
                    double l = m.get<double>();
                    setlen( {l,l,l} );
                }
                if (m.is_array()) {
                    len << m.get<vector<double>>();
                    setlen( len );
                }
                if (getVolume()>0)
                    return;
            }
        }
        catch(std::exception& e) {
            throw std::runtime_error(name + ": " + e.what());
        }
    }

    Cuboidslit::Cuboidslit() { name += " (XY-periodicity)"; }

    Cuboidslit::Cuboidslit( Tmjson &j ) : Cuboid(j) { name += " (XY-periodicity)"; }

    /**
     * @param length Length of the Cylinder (angstrom)
     * @param radius Radius of the Cylinder (angstrom)
     */
    Cylinder::Cylinder( double length, double radius ) : Geometrybase("Cylinder")
    {
        init(length, radius);
    }

    /**
     * The json object is scanned for the following parameters:
     *
     * Key      | Description
     * :------- | :-------------------------
     * `length` | Cylinder length [angstrom]
     * `radius` | Cylinder radius [angstrom] 
     */
    Cylinder::Cylinder( Tmjson &j ) : Geometrybase("Cylinder")
    {
        try {
            init( j.at("length"),
                  j.at("radius") );
        }
        catch(std::exception& e) {
            throw std::runtime_error(name + ": " + e.what());
        }
    }

    void Cylinder::init( double length, double radius )
    {
        name = "Cylinder (hard ends)";
        _len = length;
        setVolume(pc::pi * radius * radius * _len);
    }

    /**
     * Dummy function not be used other than for compatibility
     * reasons. Sets length to x component of vector.
     */
    void Cylinder::setlen( const Point &l )
    {
        init(l.x(), r);
    }

    void Cylinder::_setVolume( double newV )
    {
        r2 = newV / (pc::pi * _len);
        r = sqrt(r2);
        diameter = 2 * r;
        _halflen = _len / 2;
    }

    /**
     * - 3D: \f$ l \pi r^2 \f$ (default)
     * - 2D: \f$ \pi r^2 \f$
     * - 1D: \f$ l \f$
     */
    double Cylinder::_getVolume(int dim) const
    {
        if (dim==1)
            return _len;
        if (dim==2)
            return pc::pi * r2;
        return r2 * pc::pi * _len;
    }

    void Cylinder::boundary( Point &p ) const {}

    void Cylinder::randompos( Point &m )
    {
        double l = r2 + 1;
        m.z() = slump.half() * _len;
        while ( l > r2 )
        {
            m.x() = slump.half() * diameter;
            m.y() = slump.half() * diameter;
            l = m.x() * m.x() + m.y() * m.y();
        }
    }

    bool Cylinder::collision( const Point &a, double radius, collisiontype type ) const
    {
        assert((_halflen - _len / 2) < 1e-9 && "Cylinder length initialization problems");
        if ( a.z() < -_halflen )
            return true;
        if ( a.z() > _halflen )
            return true;
        if ( a.x() * a.x() + a.y() * a.y() > r2 )
            return true;
        return false;
    }

    string Cylinder::_info( char w )
    {
        std::ostringstream o;
        o << pad(SUB, w, "Length") << _halflen * 2 << textio::_angstrom << endl
          << pad(SUB, w, "Radius") << sqrt(r2) << textio::_angstrom << endl;
        return o.str();
    }

    Cuboid Cylinder::inscribe() const
    {
        Cuboid c;
        c.setlen({diameter, diameter, _len});
        return c;
    }

    PeriodicCylinder::PeriodicCylinder(
        double length, double radius ) : Cylinder(length, radius)
    {
        name = "Periodic " + name;
    }

    PeriodicCylinder::PeriodicCylinder( Tmjson &j ) : Cylinder(j)
    {
        name = "Periodic " + name;
    }

    void PeriodicCylinder::boundary( Point &a ) const
    {
        if ( std::abs(a.z()) > _halflen )
            a.z() -= _len * anint(a.z() / _len);
    }

    void QuaternionRotate::setAxis( Geometrybase &g, const Point &beg, const Point &end, double angle )
    {
        geoPtr = &g;
        origin = beg;
        angle_ = angle;
        Point u(end - beg); //Point u(end-beg);
        assert(u.squaredNorm() > 0 && "Rotation vector has zero length");
        g.boundary(u);
        u.normalize(); // make unit vector
        q = Eigen::AngleAxisd(angle, u);

        rot_mat << 0, -u.z(), u.y(), u.z(), 0, -u.x(), -u.y(), u.x(), 0;
        rot_mat =
            Eigen::Matrix3d::Identity() + rot_mat * std::sin(angle) + rot_mat * rot_mat * (1 - std::cos(angle));
    }

    FindSpace::FindSpace()
    {
        dir = Point(1, 1, 1);
        allowContainerOverlap = false;
        allowMatterOverlap = false;
    }

    CuboidNoPBC::CuboidNoPBC( ) { name += " (No PBC)"; }

    CuboidNoPBC::CuboidNoPBC( Tmjson &j ) : Cuboid(j) { name += " (No PBC)"; }
  }//namespace geometry
}//namespace
