#include "geometry.h"

Faunus::Geometry::GeometryBase::GeometryBase() {
    using namespace std::placeholders;
    boundaryFunc = nullptr;//std::bind( &GeometryBase::boundary, this, _1);
    distanceFunc = nullptr;//std::bind( &GeometryBase::vdist, this, _1, _2);
}

Faunus::Geometry::GeometryBase &Faunus::Geometry::GeometryBase::operator=(const Faunus::Geometry::GeometryBase &other) {
    if (this!=&other)
        name = other.name;
    return *this;
}

Faunus::Geometry::GeometryBase::GeometryBase(const Faunus::Geometry::GeometryBase &other) {
    if (this!=&other)
        name = other.name;
}

Faunus::Geometry::GeometryBase::~GeometryBase() {}

void Faunus::Geometry::Box::setLength(const Faunus::Point &l) {
    len = l;
    len_half = l*0.5;
    len_inv = l.cwiseInverse();
    c1 = l.y()/l.x();
    c2 = l.z()/l.x();
}

Faunus::Point Faunus::Geometry::Box::setVolume(double V, Faunus::Geometry::VolumeMethod method) {
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
            alpha = sqrt( V / pow(getVolume(), 1./3.) );
            s = { alpha, alpha, 1/(alpha*alpha) };
            setLength( len.cwiseProduct(s) );
            return s;
        default:
            throw std::runtime_error("unknown volume scaling method");
    }
    assert( fabs(getVolume()-V)<1e-6 );
    return s.cwiseQuotient(len); // this will scale any point to new volume
}

double Faunus::Geometry::Box::getVolume(int dim) const {
    assert(dim==3);
    return len.x() * len.y() * len.z();
}

const Faunus::Point &Faunus::Geometry::Box::getLength() const { return len; }

bool Faunus::Geometry::Box::collision(const Faunus::Point &a, double r) const {
    if ( std::fabs(a.z()+r) > len_half.z())
        return true;
    if ( std::fabs(a.x()+r) > len_half.x())
        return true;
    if ( std::fabs(a.y()+r) > len_half.y())
        return true;
    return false;
}

void Faunus::Geometry::Box::randompos(Faunus::Point &m, Faunus::Random &rand) const {
    m.x() = (rand()-0.5) * this->len.x();
    m.y() = (rand()-0.5) * this->len.y();
    m.z() = (rand()-0.5) * this->len.z();
}

void Faunus::Geometry::to_json(Faunus::json &j, const Faunus::Geometry::Box &b) {
    j["length"] = b.getLength();
}

void Faunus::Geometry::from_json(const Faunus::json &j, Faunus::Geometry::Box &b) {
    try {
        if (j.is_object()) {
            auto m = j.at("length");
            if (m.is_number()) {
                double l = m.get<double>();
                b.setLength( {l,l,l} );
            }
            else if (m.is_array())
                if (m.size()==3) {
                    Point len = m.get<Point>();
                    b.setLength( len );
                }
        }
        if (b.getVolume()<=0)
            throw std::runtime_error("box volume is zero or less");
    }
    catch(std::exception& e) {
        throw std::runtime_error( e.what() );
    }
}

void Faunus::Geometry::from_json(const Faunus::json &j, Faunus::Geometry::Cylinder &cyl) {
    cyl.set( j.at("radius"), j.at("length") );
}

void Faunus::Geometry::from_json(const Faunus::json &j, Faunus::Geometry::Sphere &g) {
    g.set( j.at("radius") );
}

void Faunus::Geometry::to_json(Faunus::json &j, const Faunus::Geometry::Sphere &g) {
    j["radius"] = g.getRadius();
}

void Faunus::Geometry::Cylinder::set(double radius, double length) {
    len = length;
    r = radius;
    r2 = r*r;
    diameter = 2*r;
    Box::setLength( { diameter, diameter, len } );
}

Faunus::Point Faunus::Geometry::Cylinder::setVolume(double V, Faunus::Geometry::VolumeMethod method) {
    assert(method==ISOTROPIC);
    double rold = r;
    r = std::sqrt( V / (pc::pi * len) );
    set( r, len );
    return Point(r/rold, r/rold, 1);
}

double Faunus::Geometry::Cylinder::getVolume(int dim) const {
    if (dim==1)
        return len;
    if (dim==2)
        return pc::pi * r2;
    return r2 * pc::pi * len;
}

void Faunus::Geometry::Cylinder::randompos(Faunus::Point &m, Faunus::Random &rand) const {
    double l = r2 + 1;
    m.z() = (rand()-0.5) * len;
    while ( l > r2 )
    {
        m.x() = (rand()-0.5) * diameter;
        m.y() = (rand()-0.5) * diameter;
        l = m.x() * m.x() + m.y() * m.y();
    }
}

bool Faunus::Geometry::Cylinder::collision(const Faunus::Point &a, double r) const {
    if ( a.z() < -0.5*len )
        return true;
    if ( a.z() > 0.5*len )
        return true;
    if ( a.x()*a.x() + a.y()*a.y() > r2 )
        return true;
    return false;
}

void Faunus::Geometry::Sphere::set(double radius) {
    r = radius;
    r2 = r*r;
    Box::setLength( {2*r,2*r,2*r} );
}

double Faunus::Geometry::Sphere::getRadius() const { return r; }

Faunus::Point Faunus::Geometry::Sphere::setVolume(double V, Faunus::Geometry::VolumeMethod method) {
    assert(method==ISOTROPIC);
    double rold = r;
    double r = std::cbrt(3*V/(4*pc::pi));
    set(r);
    return Point().setConstant(r/rold);
}

double Faunus::Geometry::Sphere::getVolume(int dim) const {
    if (dim==1)
        return 2*r;
    if (dim==2)
        return pc::pi * r2;
    return 4*pc::pi*r*r2/3;
}

void Faunus::Geometry::Sphere::randompos(Faunus::Point &a, Faunus::Random &rand) const {
    double diameter = 2*r;
    do {
        a.x() = (rand() - 0.5) * diameter;
        a.y() = (rand() - 0.5) * diameter;
        a.z() = (rand() - 0.5) * diameter;
    } while ( a.squaredNorm() > r2 );
}

bool Faunus::Geometry::Sphere::collision(const Faunus::Point &a, double r) const {
    return (a.squaredNorm()>r2) ? true : false;
}
