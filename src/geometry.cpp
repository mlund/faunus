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

Faunus::Point Faunus::Geometry::Box::getLength() const { return len; }

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

void Faunus::Geometry::from_json(const Faunus::json &j, Faunus::Geometry::Chameleon &g) {
    g.from_json(j);
}

void Faunus::Geometry::to_json(Faunus::json &j, const Faunus::Geometry::Chameleon &g) {
    g.to_json(j);
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

void Faunus::Geometry::Chameleon::setLength(const Faunus::Point &l) {
    len = l;
    len_half = l*0.5;
    len_inv = l.cwiseInverse();
    c1 = l.y()/l.x();
    c2 = l.z()/l.x();
}

double Faunus::Geometry::Chameleon::getVolume(int dim) const {
    switch (type) {
        case SPHERE:   return 4*pc::pi/3*radius*radius*radius;
        case CUBOID:   return len.x()*len.y()*len.z();
        case SLIT:     return len.x()*len.y()*len.z()*0.5;
        case CYLINDER: return pc::pi*radius*radius*len.z();
    }
}

Faunus::Point Faunus::Geometry::Chameleon::setVolume(double V, Faunus::Geometry::VolumeMethod method) {
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

    if (type==CYLINDER and method==ISOTROPIC) {
        double rold = radius;
        radius = std::sqrt( V / (pc::pi * len.z()) );
        setLength( { 4*radius, 4*radius, len.z() } );
        assert( fabs(getVolume()-V)<1e-6 );
        return Point(radius/rold, radius/rold, 1);
    }

    if (type==SPHERE and method==ISOTROPIC) {
        double rold = radius;
        radius = std::cbrt(3*V/(4*pc::pi));
        setLength( 4*Point(radius,radius,radius) ); // disable min. image in xyz
        assert( fabs(getVolume()-V)<1e-6 );
        return Point().setConstant(radius/rold);
    }

    throw std::runtime_error("unsupported volume move for geometry");

    return {1,1,1};
}

Faunus::Point Faunus::Geometry::Chameleon::getLength() const {
    switch (type) {
        case CUBOID:   return len;
        case SLIT:     return {len.x(), len.y(), 0.5*len.z()};
        case SPHERE:   return {2*radius,2*radius,2*radius};
        case CYLINDER: return {2*radius,2*radius,len.z()};
    }
}

void Faunus::Geometry::Chameleon::randompos(Faunus::Point &m, Faunus::Random &rand) const {
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
                d = m.x()*m.x() + m.y()*m.y();
            } while ( d>r2 );
            break;
        default:
            m.x() = (rand()-0.5) * len.x();
            m.y() = (rand()-0.5) * len.y();
            m.z() = (rand()-0.5) * len.z();
            break;
    }
}

bool Faunus::Geometry::Chameleon::collision(const Faunus::Point &a, double r) const {
    assert(std::fabs(r)<1e-6 && "collision test w. radius unimplemented.");
    double r2 = radius*radius;
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
            if ( std::fabs(a.z()) > 0.5*len_half.z()) return true;
            break;
        case CUBOID:
            if ( std::fabs(a.x()) > len_half.x()) return true;
            if ( std::fabs(a.y()) > len_half.y()) return true;
            if ( std::fabs(a.z()) > len_half.z()) return true;
            break;
    }
    return false;
}

void Faunus::Geometry::Chameleon::from_json(const Faunus::json &j) {
    len.setZero();
    auto it = names.find( j.at("type").get<std::string>() );
    if (it!=names.end()) {
        name = it->first;
        type = it->second;

        if ( type == CUBOID or type == SLIT ) {
            auto m = j.at("length");
            if (m.is_number()) {
                double l = m.get<double>();
                setLength( {l,l,l} );
            } else if (m.is_array())
                if (m.size()==3)
                    setLength( m.get<Point>() );
        }

        if ( type == SPHERE or type == CYLINDER )
            radius = j.at("radius").get<double>();

        if (type == SLIT)
            setLength( {len.x(), len.y(), 2*len.z() } );  // disable min. image in z

        if (type == SPHERE)
            setLength( 4*Point(radius,radius,radius) ); // disable min. image in xyz

        if (type == CYLINDER) {
            len.z() = j.at("length").get<double>();
            setLength( {4*radius,4*radius,len.z()} ); // disable min. image in xy
        }
    }
}

void Faunus::Geometry::Chameleon::to_json(Faunus::json &j) const {
    assert(!name.empty());
    switch (type) {
        case SPHERE:
            j = {{"radius",radius}};
            break;
        case CYLINDER:
            j = {{"radius",radius}, {"length",len.z()}};
            break;
        case SLIT:
            j = {{"length",Point(len.x(), len.y(), 0.5*len.z())}};
            break;
        case CUBOID:
            j = {{"length",len}};
            break;
    }
    j["type"] = name;
}
