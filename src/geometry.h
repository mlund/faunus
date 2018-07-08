#pragma once

#include "core.h"

/** @brief Faunus main namespace */
namespace Faunus {

    /** @brief Simulation geometries and related operations */
    namespace Geometry {

        typedef std::function<void(Point&)> BoundaryFunction;
        typedef std::function<Point(const Point&, const Point&)> DistanceFunction;

        enum VolumeMethod {ISOTROPIC, ISOCHORIC, XY, Z};

        struct GeometryBase {
            virtual Point setVolume(double, VolumeMethod=ISOTROPIC)=0; //!< Set volume
            virtual double getVolume(int=3) const=0; //!< Get volume
            virtual void randompos( Point&, Random& ) const=0; //!< Generate random position
            virtual Point vdist( const Point&, const Point& ) const=0; //!< (Minimum) distance between two points
            virtual void boundary( Point& ) const=0; //!< Apply boundary conditions
            virtual bool collision( const Point&, double=0) const=0; //!< Check if position lies within
            virtual const Point& getLength() const=0; //!< Side lengths

            inline double sqdist( const Point &a, const Point &b ) const {
                return vdist(a,b).squaredNorm();
            } //!< Squared (minimum) distance between two points

            BoundaryFunction boundaryFunc; //!< Functor for boundary()
            DistanceFunction distanceFunc; //!< Functor for vdist()

            std::string name;

            inline GeometryBase& operator=(const GeometryBase &other) {
                if (this!=&other)
                    name = other.name;
                return *this;
            } //!< Required since we *do not* want to overwrite functors (distance, boundary)

            inline GeometryBase(const GeometryBase &other) {
                if (this!=&other)
                    name = other.name;
            } //!< Required since we *do not* want to overwrite functors (distance, boundary)

            GeometryBase() {
                using namespace std::placeholders;
                boundaryFunc = nullptr;//std::bind( &GeometryBase::boundary, this, _1);
                distanceFunc = nullptr;//std::bind( &GeometryBase::vdist, this, _1, _2);
            }

        }; //!< Base class for all geometries

        /** @brief Cuboidal box */
        class Box : public GeometryBase {
            protected:
                Point len, //!< side length
                      len_half, //!< half side length
                      len_inv; //!< inverse side length
                double c1, c2;

            public:

                void setLength( const Point &l ) {
                    len = l;
                    len_half = l*0.5;
                    len_inv = l.cwiseInverse();
                    c1 = l.y()/l.x();
                    c2 = l.z()/l.x();
                } //!< Set cuboid side length

                Point setVolume(double V, VolumeMethod method=ISOTROPIC) override {
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
                    assert( fabs(getVolume()-V)<1e-9 );
                    return s.cwiseQuotient(len); // this will scale any point to new volume
                }

                double getVolume(int dim=3) const override {
                    assert(dim==3);
                    return len.x() * len.y() * len.z();
                }

                const Point& getLength() const override { return len; } //!< Side lengths

                bool collision(const Point &a, double r=0) const override {
                    if ( std::fabs(a.z()+r) > len_half.z())
                    return true;
                    if ( std::fabs(a.x()+r) > len_half.x())
                    return true;
                    if ( std::fabs(a.y()+r) > len_half.y())
                    return true;
                    return false;
                }

                void randompos( Point &m, Random &rand ) const override {
                    m.x() = (rand()-0.5) * this->len.x();
                    m.y() = (rand()-0.5) * this->len.y();
                    m.z() = (rand()-0.5) * this->len.z();
                }
        };

        void to_json(json &j, const Box &b) {
            j["length"] = b.getLength();
        }

        void from_json(const json& j, Box &b) {
            try {
                if (j.is_object()) {
                    auto m = j.at("length");
                    if (m.is_number()) {
                        double l = m.get<double>();
                        b.setLength( {l,l,l} );
                    }
                    else if (m.is_array()) {
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

        /** @brief Periodic boundary conditions */
        template<bool X=true, bool Y=true, bool Z=true>
            struct PBC : public Box {

                PBC() {
                    using namespace std::placeholders;
                    boundaryFunc = std::bind( &PBC<X,Y,Z>::boundary, this, _1);
                    distanceFunc = std::bind( &PBC<X,Y,Z>::vdist, this, _1, _2);
                    distanceFunc = [this](const Point &i, const Point &j){return this->vdist(i,j);};
                    boundaryFunc = [this](Point &i){this->boundary(i);};
                }

                Point vdist( const Point &a, const Point &b ) const override
                {
                    Point r(a - b);
                    if (X) {
                        if ( r.x() > len_half.x())
                            r.x() -= len.x();
                        else if ( r.x() < -len_half.x())
                            r.x() += len.x();
                    }
                    if (Y) {
                        if ( r.y() > len_half.y())
                            r.y() -= len.y();
                        else if ( r.y() < -len_half.y())
                            r.y() += len.y();
                    }
                    if (Z) {
                        if ( r.z() > len_half.z())
                            r.z() -= len.z();
                        else if ( r.z() < -len_half.z())
                            r.z() += len.z();
                    }
                    return r;
                } //!< (Minimum) distance between two points

                void unwrap( Point &a, const Point &ref ) const {
                    a = vdist(a, ref) + ref;
                } //!< Remove PBC with respect to a reference point

                template<typename T=double>
                    inline int anint( T x ) const
                    {
                        return int(x > 0.0 ? x + 0.5 : x - 0.5);
                    } //!< Round to int

                void boundary( Point &a ) const override
                {
                    if (X)
                        if ( std::fabs(a.x()) > len_half.x())
                            a.x() -= len.x() * anint(a.x() * len_inv.x());
                    if (Y)
                        if ( std::fabs(a.y()) > len_half.y())
                            a.y() -= len.y() * anint(a.y() * len_inv.y());
                    if (Z)
                        if ( std::fabs(a.z()) > len_half.z())
                            a.z() -= len.z() * anint(a.z() * len_inv.z());
                } //!< Apply boundary conditions

            };

        using Cuboid = PBC<true,true,true>; //!< Cuboid w. PBC in all directions
        using Cuboidslit = PBC<true,true,false>; //!< Cuboidal slit w. PBC in XY directions
        using CuboidNoPBC = PBC<false,false,false>; //!< Hard cuboid - no PBC

#ifdef DOCTEST_LIBRARY_INCLUDED
        TEST_CASE("[Faunus] PBC/Cuboid") {
            Cuboid geo = R"( {"length": [2,3,4]} )"_json;

            CHECK( geo.getVolume() == doctest::Approx(2*3*4) );

            Point a(1.1, 1.5, -2.001);
            geo.boundaryFunc(a);
            CHECK( a.x() == doctest::Approx(-0.9) );
            CHECK( a.y() == doctest::Approx(1.5) );
            CHECK( a.z() == doctest::Approx(1.999) );
            Point b = a;
            geo.boundary(b);
            CHECK( a == b );

            // check copying, in particular functors
            auto g1 = new Cuboid();
            g1->boundaryFunc(a);
            geo = *g1; // functors should *not be copied*
            delete g1;
            geo.boundaryFunc(a);
        }
#endif

        /** @brief Cylindrical cell */
        class Cylinder : public PBC<false,false,true> {
            private:
                double r, r2, diameter, len;
                typedef PBC<false,false,true> base;
            public:
                void set(double radius, double length) {
                    len = length;
                    r = radius;
                    r2 = r*r;
                    diameter = 2*r;
                    Box::setLength( { diameter, diameter, len } );
                } //!< Set radius

                Point setVolume(double V, VolumeMethod method=ISOTROPIC) override {
                    assert(method==ISOTROPIC);
                    double rold = r;
                    r = std::sqrt( V / (pc::pi * len) );
                    set( r, len );
                    return Point(r/rold, r/rold, 1);
                } //!< Adjust radius to match volume; length is conserved.

                double getVolume(int dim=3) const override {
                    if (dim==1)
                        return len;
                    if (dim==2)
                        return pc::pi * r2;
                    return r2 * pc::pi * len;
                } //<! Return volume

                void randompos( Point &m, Random &rand ) const override
                {
                    double l = r2 + 1;
                    m.z() = (rand()-0.5) * len;
                    while ( l > r2 )
                    {
                        m.x() = (rand()-0.5) * diameter;
                        m.y() = (rand()-0.5) * diameter;
                        l = m.x() * m.x() + m.y() * m.y();
                    }
                }

                bool collision(const Point &a, double r=0) const override {
                    if ( a.z() < -0.5*len )
                        return true;
                    if ( a.z() > 0.5*len )
                        return true;
                    if ( a.x()*a.x() + a.y()*a.y() > r2 )
                        return true;
                    return false;
                }
        };

        void from_json(const json& j, Cylinder &cyl) {
            cyl.set( j.at("radius"), j.at("length") );
        }

#ifdef DOCTEST_LIBRARY_INCLUDED
        TEST_CASE("[Faunus] Cylinder") {
            Cylinder c;
            c.set( 1.0, 1/pc::pi );
            CHECK( c.getVolume() == doctest::Approx( 1.0 ) );
            CHECK( c.collision({1.01,0,0}) == true);
            CHECK( c.collision({0.99,0,0}) == false);
            CHECK( c.collision({0.99,1/pc::pi+0.01,0}) == true);
        }
#endif

        /** @brief Spherical cell */
        class Sphere : public PBC<false,false,false> {
            private:
                double r, r2;
            public:
                void set(double radius) {
                    r = radius;
                    r2 = r*r;
                    Box::setLength( {2*r,2*r,2*r} );
                } //!< Set radius

                double getRadius() const { return r; }

                Point setVolume(double V, VolumeMethod method=ISOTROPIC) override {
                    assert(method==ISOTROPIC);
                    double rold = r;
                    double r = std::cbrt(3*V/(4*pc::pi));
                    set(r);
                    return Point().setConstant(r/rold);
                } //!< Adjust radius to match volume; length is conserved.

                double getVolume(int dim=3) const override {
                    if (dim==1)
                        return 2*r;
                    if (dim==2)
                        return pc::pi * r2;
                    return 4*pc::pi*r*r2/3;
                } //<! Return volume

                void randompos( Point &a, Random &rand ) const override {
                    double diameter = 2*r;
                    do {
                        a.x() = (rand() - 0.5) * diameter;
                        a.y() = (rand() - 0.5) * diameter;
                        a.z() = (rand() - 0.5) * diameter;
                    } while ( a.squaredNorm() > r2 );
                }

                bool collision(const Point &a, double r=0) const override {
                    return (a.squaredNorm()>r2) ? true : false;
                }
        };

        void from_json(const json& j, Sphere &g) {
            g.set( j.at("radius") );
        }

        void to_json(json &j, const Sphere &g) {
            j["radius"] = g.getRadius();
        }

#ifdef DOCTEST_LIBRARY_INCLUDED
        TEST_CASE("[Faunus] Sphere") {
            Sphere c = R"( { "radius" : 2 } )"_json;
            CHECK( c.getVolume() == doctest::Approx( 4*pc::pi*8/3.0 ) );
            CHECK( c.collision( { 2.01, 0, 0 } ) == true );
            CHECK( c.collision( { 1.99, 0, 0 } ) == false );
            CHECK( c.getLength().squaredNorm() == doctest::Approx(48) );
            c.setVolume(123.4);
            CHECK( c.getVolume() == doctest::Approx(123.4) );
        }
#endif

        enum class weight { MASS, CHARGE, GEOMETRIC };

        template<typename Titer, typename Tparticle=typename Titer::value_type, typename weightFunc>
            Point anyCenter( Titer begin, Titer end, BoundaryFunction boundary, const weightFunc &weight,
                    const Point &shift={0,0,0})
            {
                double sum=0;
                Point c(0,0,0);
                try {
                    for (auto &i=begin; i!=end; ++i) {
                        Point t = i->pos + shift;       // translate to origo
                        boundary(t);
                        double w = weight(*i);
                        c += w * t;
                        sum += w;
                    }
                    c = c/sum - shift;
                    boundary(c);
                    if (std::isnan(c[0]))
                        throw std::runtime_error("sum of weights is zero");
                }
                catch(std::exception& e) {
                    throw std::runtime_error("anyCenter error: " + std::string(e.what()));
                }
                return c;
            } //!< Mass, charge, or geometric center of a collection of particles

        template<typename Titer, typename Tparticle=typename Titer::value_type>
            Point massCenter(Titer begin, Titer end, BoundaryFunction boundary=[](Point&){}, const Point &shift={0,0,0}) {
                return anyCenter(begin, end, boundary,
                        []( const Tparticle &p ){ return atoms<Tparticle>.at(p.id).mw; }, shift );
            } // Mass center

#ifdef DOCTEST_LIBRARY_INCLUDED
        TEST_CASE("[Faunus] anyCenter") {
            typedef Particle<Radius, Charge, Dipole, Cigar> T;

            Cylinder cyl = json( {{"length",100}, {"radius",20}} );
            std::vector<T> p;

            CHECK( !atoms<T>.empty() ); // set in a previous test
            p.push_back( atoms<T>[0].p );
            p.push_back( atoms<T>[0].p );

            p.front().pos = {10, 10, -10};
            p.back().pos  = {15, -10, 10};

            Point cm = Geometry::massCenter(p.begin(), p.end(), cyl.boundaryFunc);
            CHECK( cm.x() == doctest::Approx(12.5) );
            CHECK( cm.y() == doctest::Approx(0) );
            CHECK( cm.z() == doctest::Approx(0) );
        }
#endif

        template<class Titer=typename std::vector<T>::iterator>
            void translate( Titer begin, Titer end, const Point &d,
                    BoundaryFunction boundary=[](Point &i){}  )
            {
                for ( auto i=begin; i!=end; ++i )
                {
                    i->pos += d;
                    boundary(i->pos);
                }
            } //!< Vector displacement of a range of particles

        template<typename Titer>
            void cm2origo( Titer begin, Titer end, BoundaryFunction boundary=[](Point&){} )
            {
                Point cm = massCenter(begin, end, boundary);
                translate(begin, end, -cm, boundary);
            } //!< Translate to that mass center is in (0,0,0)

        template<typename Titer>
            void rotate(
                    Titer begin,
                    Titer end,
                    const Eigen::Quaterniond &q,
                    BoundaryFunction boundary=[](Point&){},
                    const Point &shift=Point(0,0,0) )
            {
                auto m = q.toRotationMatrix(); // rotation matrix
                for (auto i=begin; i!=end; ++i) {
                    i->rotate(q,m); // rotate internal coordinates
                    i->pos += shift;
                    boundary(i->pos);
                    i->pos = q*i->pos;
                    boundary(i->pos);
                    i->pos -= shift;
                    boundary(i->pos);
                }
            } //!< Rotate particle pos and internal coordinates

        /* 
         * @brief Calculate mass center of cluster of particles in unbounded environment 
         *
         * [More info](http://dx.doi.org/10.1080/2151237X.2008.10129266)
         */
        template<class Tspace, class GroupIndex>
            Point trigoCom( const Tspace &spc, const GroupIndex &groups, const std::vector<int> &dir = {0, 1, 2} )
            {
                assert(!dir.empty() && dir.size() <= 3);
                Point xhi(0, 0, 0), zeta(0, 0, 0), theta(0, 0, 0), com(0, 0, 0);
                for ( auto k : dir ) {
                    double q = 2 * pc::pi / spc.geo.getLength()[k];
                    size_t N=0;
                    for (auto i : groups)
                        for (auto &particle : spc.groups[i]) {
                            theta[k] = particle.pos[k] * q;
                            zeta[k] += std::sin(theta[k]);
                            xhi[k] += std::cos(theta[k]);
                            N++;
                        }
                    theta[k] = std::atan2(-zeta[k] / N, -xhi[k] / N) + pc::pi;
                    com[k] = spc.geo.getLength()[k] * theta[k] / (2 * pc::pi);
                }
                spc.geo.boundary(com); // is this really needed?
                return com;
            }

        /**
         * @brief Calculates the gyration tensor of a molecular group
         * The gyration tensor is computed from the dyadic product of the position
         * vectors in the c-o-m reference system, \f$ t_{i} = r_{i} - shift \f$:
         * \f$ S = \sum_{i=0}^{N} t_{i} t_{i}^{T} \f$
         */
        template<typename iter>
            Tensor gyration(iter begin, iter end, BoundaryFunction boundary=[](const Point&){}, const Point shift=Point(0,0,0) ) {
                Tensor S;
                size_t n = std::distance(begin,end);
                if (n>0) {
                    for (auto it=begin; it!=end; ++it) {
                        Point t = it->pos - shift;
                        boundary(t);
                        // S += t * t.transpose();
                        S += t.squaredNorm() * Eigen::Matrix<double, 3, 3>::Identity() - t * t.transpose();
                    }
                    return S*(1.0/n);
                }
                return S;
            }

        template<class Titer>
            double monopoleMoment( Titer begin, Titer end ) {
                double z=0;
                for (auto it=begin; it!=end; ++it)
                    z += it->charge;
                return z;
            } //!< Calculates dipole moment vector

        template<class Titer>
            Point dipoleMoment( Titer begin, Titer end, BoundaryFunction boundary=[](const Point&){}, double cutoff=pc::infty) {
                Point mu(0,0,0);
                for (auto it=begin; it!=end; ++it) {
                    Point t = it->pos - begin->pos;
                    boundary(t);
                    if (t.squaredNorm() < cutoff*cutoff)
                        mu += t * it->charge;
                }
                return mu;
            } //!< Calculates dipole moment vector

        template<class Titer>
            Tensor quadrupoleMoment( Titer begin, Titer end, BoundaryFunction boundary=[](const Point&){}, double cutoff=pc::infty) {
                Tensor theta;
                theta.setZero();
                for (auto it=begin; it!=end; ++it) {
                    Point t = it->pos - begin->pos;
                    boundary(t);
                    if (t.squaredNorm() < cutoff*cutoff)
                        theta += t * t.transpose() * it->charge;
                }
                return 0.5 * theta;
            } //!< Calculates quadrupole moment tensor (with trace)

        template<class Tgroup>
            auto toMultipole(const Tgroup &g, BoundaryFunction boundary=[](const Point&){}, double cutoff=pc::infty) {
                Particle<Charge,Dipole,Quadrupole> m;
                m.pos = g.cm;
                m.charge = Geometry::monopoleMoment(g.begin(), g.end());                   // monopole
                m.mu = Geometry::dipoleMoment(g.begin(), g.end(), boundary, cutoff);   // dipole
                m.Q = Geometry::quadrupoleMoment(g.begin(), g.end(), boundary, cutoff);// quadrupole
                m.mulen = m.mu.norm();
                if (m.mulen>1e-9)
                    m.mu.normalize();
                return m;
            } //<! Group --> Multipole

    } //geometry namespace
}//end of faunus namespace
