#include "geometry.h"

namespace Faunus {
    Point ranunit(Random &rand, const Point &dir) {
        double r2;
        Point p;
        do {
            for (size_t i=0; i<3; i++)
                p[i] = ( rand()-0.5 ) * dir[i];
            r2 = p.squaredNorm();
        } while ( r2 > 0.25 );
        return p / std::sqrt(r2);
    }

    Point ranunit_polar(Random &rand) {
        return rtp2xyz( {1, 2*pc::pi*rand(), std::acos(2*rand()-1)} );
    }

    Point xyz2rth(const Point &p, const Point &origin, const Point &dir, const Point &dir2) {
        assert( fabs(dir.norm()-1.0)<1e-6 );
        assert( fabs(dir2.norm()-1.0)<1e-6 );
        assert( fabs(dir.dot(dir2))<1e-6 ); // check if unit-vectors are perpendicular
        Point xyz = p - origin;
        double h = xyz.dot(dir);
        Point xy = xyz - dir*h;
        double x = xy.dot(dir2);
        Point y = xy - dir2*x;
        double theta = std::atan2(y.norm(),x);
        double radius = xy.norm();
        return {radius,theta,h};
    }

    Point xyz2rtp(const Point &p, const Point &origin) {
        Point xyz = p - origin;
        double radius = xyz.norm();
        return {
            radius,
                std::atan2( xyz.y(), xyz.x() ),
                std::acos( xyz.z()/radius) };
    }

    Point rtp2xyz(const Point &rtp, const Point &origin) {
        return origin + rtp.x() * Point(
                std::cos(rtp.y()) * std::sin(rtp.z()),
                std::sin(rtp.y()) * std::sin(rtp.z()),
                std::cos(rtp.z()) );
    }

    namespace Geometry {

        // =============== GeometryBase ===============

        GeometryBase::~GeometryBase() {}

        const std::map<std::string, Variant> GeometryBase::names = {
                {
                        {"cuboid", CUBOID},
                        {"cylinder", CYLINDER},
                        {"slit", SLIT},
                        {"sphere", SPHERE},
                        {"hexagonal", HEXAGONAL},
                        {"octahedron", OCTAHEDRON}
                }
        };


        // =============== Cuboid ===============

        Cuboid::Cuboid(const Point &p) {
            setLength(p);
        }

        Cuboid::Cuboid(double x, double y, double z) : Cuboid(Point(x, y, z)) {
        }

        Cuboid::Cuboid(double x) : Cuboid(x, x, x) {
        }

        inline Point Cuboid::getLength() const {
            return box;
        }

        void Cuboid::setLength(const Point &len) {
            box = len;
            box_half = 0.5 * box;
            box_inv = box.cwiseInverse();
            ratio_yx = box.y() / box.x();
            ratio_zx = box.z() / box.x();
        }

        inline double Cuboid::getVolume(int) const {
            return box.x() * box.y() * box.z();
        }

        Point Cuboid::setVolume(double volume, const VolumeMethod method) {
            //TODO a better implementation shall assign scaling and compute the new box as a product
            //TODO move out "isochoric" scaling
            double x, alpha;
            Point new_box, box_scaling;
            switch (method) {
                case ISOTROPIC:
                    x = std::cbrt(volume / (ratio_yx * ratio_zx)); // keep aspect ratios
                    new_box = {x, ratio_yx * x, ratio_zx * x};
                    break;
                case XY:
                    x = std::sqrt(volume / (ratio_yx * box.z())); // keep xy aspect ratio
                    new_box = {x, ratio_yx * x, box.z()};         // z is untouched
                    break;
                case ISOCHORIC:
                    // z is scaled by 1/alpha/alpha, x and y are scaled by alpha
                    alpha = std::cbrt(volume / (ratio_yx * ratio_zx)) / box.x();
                    new_box = box.cwiseProduct(Point(alpha, alpha, 1 / (alpha * alpha)));
                    break;
                default:
                    throw std::invalid_argument("unsupported volume scaling method for the cuboid geometry");
            }
            box_scaling = new_box.cwiseQuotient(box);
            setLength(new_box);
            assert(fabs(getVolume() - volume) < 1e-6);
            return box_scaling; // this will scale any point to new volume
        }

        inline void Cuboid::boundary(Point &a) const {
            // xyz-pbc
            if (std::fabs(a.x()) > box_half.x())
                a.x() -= box.x() * anint(a.x() * box_inv.x());
            if (std::fabs(a.y()) > box_half.y())
                a.y() -= box.y() * anint(a.y() * box_inv.y());
            if (std::fabs(a.z()) > box_half.z())
                a.z() -= box.z() * anint(a.z() * box_inv.z());
        }

        inline Point Cuboid::vdist(const Point &a, const Point &b) const {
            // xyz-pbc
            Point distance(a - b);
            if (distance.x() > box_half.x())
                distance.x() -= box.x();
            else if (distance.x() < -box_half.x())
                distance.x() += box.x();
            if (distance.y() > box_half.y())
                distance.y() -= box.y();
            else if (distance.y() < -box_half.y())
                distance.y() += box.y();
            if (distance.z() > box_half.z())
                distance.z() -= box.z();
            else if (distance.z() < -box_half.z())
                distance.z() += box.z();
            return distance;
        }

        inline void Cuboid::randompos(Point &m, Random &rand) const {
            m.x() = (rand() - 0.5) * box.x();
            m.y() = (rand() - 0.5) * box.y();
            m.z() = (rand() - 0.5) * box.z();
        }

        inline bool Cuboid::collision(const Point &a) const {
            bool collision =
                    std::fabs(a.x()) > box_half.x() ||
                    std::fabs(a.y()) > box_half.y() ||
                    std::fabs(a.z()) > box_half.z();
            return collision;
        }

        void Cuboid::from_json(const json &j) {
            std::tie(std::ignore, name) = variantName(j);
            box.setZero();

            auto m = j.at("length");
            if (m.is_number()) {
                double l = m.get<double>();
                setLength({l, l, l});
            } else if (m.is_array()) {
                if (m.size() == 3) {
                    setLength(m.get<Point>());
                } else {
                    // TODO warning
                }
            } else {
                // TODO warning
            }
        }

        void Cuboid::to_json(json &j) const {
            j = {{"type",   name},
                 {"length", box}};
        }


        // =============== Chameleon==============

        void from_json(const json &j, Chameleon &g) {
            try { 
                g.from_json(j);
            } catch(std::exception& e) {
                throw std::runtime_error("geometry construction error: "s + e.what() + usageTip["geometry"]);
            }
        }

        void to_json(json &j, const Chameleon &g) {
            g.to_json(j);
        }

        void Chameleon::setLength(const Point &l) {
            len = l;
            len_half = l*0.5;
            len_inv = l.cwiseInverse();
            c1 = l.y()/l.x();
            c2 = l.z()/l.x();
        }

        double Chameleon::getVolume(int) const {
            switch (type) {
                case SPHERE:     return 4*pc::pi/3*radius*radius*radius;
                case CUBOID:     return len.x()*len.y()*len.z();
                case SLIT:       return len.x()*len.y()*len.z() / pbc_disable;
                case CYLINDER:   return pc::pi*radius*radius*len.z();
                case HEXAGONAL:  return 2.0*std::sqrt(3.0)*radius*radius*len.z();
                case OCTAHEDRON: return 8.0*std::sqrt(2.0)*radius*radius*radius; // for the truncated octahedron then 'radius' really is the side-length 'a'
            }
            assert(false);
            return 0;
        }

        Point Chameleon::setVolume(double V, VolumeMethod method) {
            if (type==CUBOID) {
                double x, alpha;
                Point s;
                switch (method) {
                    case ISOTROPIC:
                        x = std::cbrt( V / (c1*c2) ); // keep aspect ratio
                        s = { x, c1*x, c2*x };
                        setLength(s);
                        break;
                    case XY:
                        x = std::sqrt( V / (c1*len.z()) ); // keep xy aspect ratio
                        s = { x, c1*x, len.z() };
                        setLength(s);  // z is untouched
                        break;
                    case ISOCHORIC:
                        // z is scaled by 1/alpha/alpha, x and y are scaled by alpha
                        alpha = std::cbrt( V / (c1*c2) ) / len.x();
                        s = { alpha, alpha, 1/(alpha*alpha) };
                        setLength( len.cwiseProduct(s) );
                        return s;
                    default:
                        throw std::runtime_error("unknown volume scaling method");
                }
                assert( fabs(getVolume()-V)<1e-6 );
                return s.cwiseQuotient(len); // this will scale any point to new volume
            }

            if (type==HEXAGONAL) {
                double alpha, oldradius, ratio;
                Point s;
                switch (method) {
                    case ISOTROPIC:
                        oldradius = radius;
                        radius = std::cbrt( V / 4.0 / sqrt(3.0) / c2 ); // keep aspect ratio, note that c2 = h / ( 2*radius )
                        setLength( { 2.0*radius, 2.0*radius, c2*2.0*radius } );
                        ratio = radius/oldradius;
                        return {ratio,ratio,ratio}; // last term comes from keeping the aspect ratio
                    case XY:
                        oldradius = radius;
                        radius = std::sqrt(V/len.z()/2.0/std::sqrt(3.0)); // inner radius
                        setLength( { 2.0*radius, 2.0*radius, len.z() } );
                        return {radius/oldradius, radius/oldradius, 1};
                    case ISOCHORIC:
                        // z is scaled by 1/alpha/alpha, radius is scaled by alpha
                        alpha = std::cbrt( V / 4.0 / sqrt(3.0) / c2 ) / radius;
                        s = { alpha, alpha, 1/(alpha*alpha) };
                        radius = alpha*radius;
                        setLength( len.cwiseProduct(s) );
                        return s;
                    default:
                        throw std::runtime_error("unknown volume scaling method");
                }
                assert( fabs(getVolume()-V)<1e-6 );
            }

            if (type==OCTAHEDRON and method==ISOTROPIC) {
                const double oldradius = radius; // for the truncated octahedron then 'radius' really is the side-length 'a'
                radius = std::cbrt( V / 8.0 / std::sqrt(2.0) );
                setLength( 2.0*radius*Point(std::sqrt(1.5),std::sqrt(2.0),std::sqrt(10.0)/2.0) ); // 2*( origin to regular hexagon, origin to square, and circumradius)
                assert( fabs(getVolume()-V)<1e-6 );
                const double ratio = radius/oldradius;
                return {ratio,ratio,ratio};
            }

            if (type==CYLINDER) {
                    double oldradius, alpha;
                    Point s;
                switch (method) {
                    case ISOTROPIC:
                        oldradius = radius;
                        radius = std::sqrt( V / (pc::pi * len.z()) );
                        setLength( { pbc_disable*2*radius, pbc_disable*2*radius, len.z() } );
                        assert( fabs(getVolume()-V)<1e-6 );
                        return {radius/oldradius, radius/oldradius, 1};
                    case ISOCHORIC:
                        // z is scaled by 1/alpha/alpha, radius is scaled by alpha
                        alpha = std::sqrt( V / (pc::pi * len.z()) ) / radius;
                        s = { alpha, alpha, 1/(alpha*alpha) };
                        radius = alpha*radius;
                        setLength( len.cwiseProduct(s) );
                        return s;
                    default:
                        throw std::runtime_error("unsupported volume move for geometry");
                }
            }

            if (type==SPHERE and method==ISOTROPIC) {
                const double oldradius = radius;
                radius = std::cbrt(3*V/(4*pc::pi));
                setLength( pbc_disable*2*Point(radius,radius,radius) ); // disable min. image in xyz
                assert( std::fabs(getVolume()-V)<1e-6 && "error setting sphere volume");
                return Point().setConstant(radius/oldradius);
            }

            throw std::runtime_error("unsupported volume move for geometry");

            return {1,1,1};
        }

        Point Chameleon::getLength() const {
            switch (type) {
                case CUBOID:     return len;
                case SLIT:       return {len.x(), len.y(), len.z() / pbc_disable};
                case SPHERE:     return {2*radius,2*radius,2*radius};
                case CYLINDER:   return {2*radius,2*radius,len.z()};
                case HEXAGONAL:  return {radius,radius,len.z()};
                case OCTAHEDRON: return {2.0*std::sqrt(1.5)*radius,2.0*std::sqrt(2.0)*radius,2.0*std::sqrt(10.0)/2.0*radius}; // 2*( origin to regular hexagon, origin to square, and circumradius )
            }
            assert(false);
            return Point();
        } //!< Enscribe geometry in cuboid

        void Chameleon::randompos(Point &m, Random &rand) const {
            double r2 = radius*radius, d=2*radius;
            switch (type) {
                case SPHERE:
                    do {
                        m.x() = (rand() - 0.5) * d;
                        m.y() = (rand() - 0.5) * d;
                        m.z() = (rand() - 0.5) * d;
                    } while ( m.squaredNorm() > r2 );
                    break;
                case CYLINDER:
                    m.z() = (rand()-0.5) * len.z();
                    do {
                        m.x() = (rand()-0.5) * d;
                        m.y() = (rand()-0.5) * d;
                    } while ( m.x()*m.x() + m.y()*m.y() > r2 );
                    break;
                case SLIT:
                    m.x() = (rand()-0.5) * len.x();
                    m.y() = (rand()-0.5) * len.y();
                    m.z() = (rand()-0.5) * len.z() / pbc_disable;
                    break;
                case CUBOID:
                    m.x() = (rand()-0.5) * len.x();
                    m.y() = (rand()-0.5) * len.y();
                    m.z() = (rand()-0.5) * len.z();
                    break;
                case OCTAHEDRON:
                    d = len.z(); // use circumradius
                    r2 = len_half.z()*len_half.z();
                    do {
                        do {
                            m.x() = (rand() - 0.5)*d;
                            m.y() = (rand() - 0.5)*d;
                            m.z() = (rand() - 0.5)*d;
                        } while ( m.squaredNorm() > r2 );
                    } while(Chameleon::collision(m));
                    break;
                case HEXAGONAL:
                    double Ra = rand();
                    double Rb = rand();
                    const double sqrtThree = sqrt(3.0);
                    const Point unitvA = Point(sqrtThree/2.0,0.5,0.0);
                    const Point unitvB = Point(sqrtThree/2.0,-0.5,0.0);
                    double R = d/sqrtThree;
                    Point p = R*unitvA*Ra+R*unitvB*Rb;
                    if(2.0*p.x() > d) p.x() = p.x() - d;
                    double theta = pc::pi/3.0*std::floor(6.0*rand());
                    double cosT = std::cos(theta);
                    double sinT = std::sin(theta);
                    m.x() = cosT*p.x() - sinT*p.y();
                    m.y() = cosT*p.y() + sinT*p.x();
                    m.z() = (rand() - 0.5)*len.z();
                    break;
            }
        }

        bool Chameleon::collision(const Point &a) const {
            double r2 = radius*radius;
            double sqrtThreeByTwo = std::sqrt(3.0)/2.0;
            double sqrtThreeI = 1.0/std::sqrt(3.0);
            switch (type) {
                case SPHERE:
                    return (a.squaredNorm()>r2) ? true : false;
                    break;
                case CYLINDER:
                    if ( std::fabs(a.z()) > len_half.z()) return true;
                    if ( a.x()*a.x() + a.y()*a.y() > r2 ) return true;
                    break;
                case SLIT:
                    if ( std::fabs(a.x()) > len_half.x()) return true;
                    if ( std::fabs(a.y()) > len_half.y()) return true;
                    if ( std::fabs(a.z()) > len_half.z() / pbc_disable ) return true;
                    break;
                case CUBOID:
                    if ( std::fabs(a.x()) > len_half.x()) return true;
                    if ( std::fabs(a.y()) > len_half.y()) return true;
                    if ( std::fabs(a.z()) > len_half.z()) return true;
                    break;
                case HEXAGONAL:
                    if(std::fabs(a.z()) > len_half.z()) return true;
                    if(std::fabs(a.dot(Point(1.0,0.0,0.0))) > len_half.x()) return true;
                    if(std::fabs(a.dot(Point(0.5,sqrtThreeByTwo,0.0))) > len_half.x()) return true;
                    if(std::fabs(a.dot(Point(-0.5,sqrtThreeByTwo,0.0))) > len_half.x()) return true;
                    break;
                case OCTAHEDRON:
                    if(std::fabs(a.dot(Point(1.0,0.0,0.0))) > len_half.y()) return true;
                    if(std::fabs(a.dot(Point(0.0,1.0,0.0))) > len_half.y()) return true;
                    if(std::fabs(a.dot(Point(0.0,0.0,1.0))) > len_half.y()) return true;
                    if(std::fabs(a.dot(Point(1.0,1.0,1.0)*sqrtThreeI)) > len_half.x()) return true;
                    if(std::fabs(a.dot(Point(1.0,1.0,-1.0)*sqrtThreeI)) > len_half.x()) return true;
                    if(std::fabs(a.dot(Point(1.0,-1.0,-1.0)*sqrtThreeI)) > len_half.x()) return true;
                    if(std::fabs(a.dot(Point(1.0,-1.0,1.0)*sqrtThreeI)) > len_half.x()) return true;
                    break;
            }
            return false;
        }

        void Chameleon::from_json(const json &j) {
            auto it = names.find( j.at("type").get<std::string>() );
            if (it==names.end())
                throw std::runtime_error("unknown geometry");

            std::tie(type, name) = variantName(j);
            len.setZero();

            if ( type == CUBOID or type == SLIT ) {
                auto m = j.at("length");
                if (m.is_number()) {
                    double l = m.get<double>();
                    setLength( {l,l,l} );
                } else if (m.is_array())
                    if (m.size()==3)
                        setLength( m.get<Point>() );
            }

            if ( type == SPHERE or type == CYLINDER or type == HEXAGONAL or type == OCTAHEDRON )
                radius = j.at("radius").get<double>();

            if (type == SLIT)
                setLength( {len.x(), len.y(), pbc_disable*len.z() } );  // disable min. image in z

            if (type == SPHERE)
                setLength( pbc_disable*2*Point(radius,radius,radius) ); // disable min. image in xyz

            if (type == OCTAHEDRON)
                setLength( {2.0*std::sqrt(1.5)*radius,2.0*std::sqrt(2.0)*radius,2.0*std::sqrt(10.0)/2.0*radius} ); // 2*( origin to regular hexagon, origin to square, and circumradius )

            if (type == CYLINDER) {
                len.z() = j.at("length").get<double>();
                setLength( {2.0*radius*pbc_disable, 2.0*radius*pbc_disable, len.z() } ); // disable min. image in xy
            }
            if (type == HEXAGONAL) {
                len.z() = j.at("length").get<double>();
                setLength( {2.0*radius, 2.0*radius, len.z() } );
            }
        }

        void Chameleon::to_json(json &j) const {
            assert(!name.empty());
            switch (type) {
                case SPHERE:
                    j = {{"radius",radius}};
                    break;
                case CYLINDER:
                    j = {{"radius",radius}, {"length",len.z()}};
                    break;
                case SLIT:
                    j = {{"length",Point(len.x(), len.y(), len.z() / pbc_disable)}};
                    break;
                case CUBOID:
                    j = {{"length",len}};
                    break;
                case HEXAGONAL:
                    j = {{"radius",radius}, {"length",len.z()}};
                    break;
                case OCTAHEDRON:
                    j = {{"radius",radius}};
                    break;
            }
            j["type"] = name;
        }

    } // end of Geometry namespace

} // end of Faunus namespace
