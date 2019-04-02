#pragma once

#include "core.h"
#include "atomdata.h"
#include "particle.h"

/** @brief Faunus main namespace */
namespace Faunus {

    /**
     * @brief Convert cartesian- to spherical-coordinates
     * @note Input (x,y,z), output \f$ (r,\theta,\phi) \f$  where \f$ r\in [0,\infty) \f$, \f$ \theta\in [-\pi,\pi) \f$, and \f$ \phi\in [0,\pi] \f$.
     */
    Point xyz2rtp(const Point&, const Point &origin={0,0,0});

    /**
     * @brief Convert spherical- to cartesian-coordinates
     * @param origin The origin to be added (optional)
     * @note Input \f$ (r,\theta,\phi) \f$  where \f$ r\in [0,\infty) \f$, \f$ \theta\in [0,2\pi) \f$, and \f$ \phi\in [0,\pi] \f$, and output (x,y,z).
     */
    Point rtp2xyz(const Point &rtp, const Point &origin = {0,0,0});

#ifdef DOCTEST_LIBRARY_INCLUDED
    TEST_CASE("[Faunus] spherical coordinates") {
        using doctest::Approx;

        Point sph1 = {2, 0.5, -0.3};
        auto pnt1 = rtp2xyz(sph1); // sph --> cart
        auto sph2 = xyz2rtp(pnt1); // cart --> sph

        CHECK( pnt1.norm() == Approx(2));
        CHECK( sph1.x() == Approx(sph2.x()));
        //CHECK( sph1.y() == Approx(sph2.y()));
        //CHECK( sph1.z() == Approx(sph2.z()));
    }
#endif

    Point ranunit(Random&, const Point &dir={1,1,1}); //!< Random unit vector using Neuman's method ("sphere picking")
    Point ranunit_polar(Random&); //!< Random unit vector using polar coordinates ("sphere picking")

#ifdef DOCTEST_LIBRARY_INCLUDED
    TEST_CASE("[Faunus] ranunit") {
        Random r;
        int n=2e5;
        Point rtp(0,0,0);
        for (int i=0; i<n; i++)
            rtp += xyz2rtp( ranunit(r) );
        rtp = rtp / n;
        CHECK( rtp.x() == doctest::Approx(1) );
        CHECK( rtp.y() == doctest::Approx(0).epsilon(0.005) ); // theta [-pi:pi] --> <theta>=0
        CHECK( rtp.z() == doctest::Approx(pc::pi/2).epsilon(0.005) );// phi [0:pi] --> <phi>=pi/2
    }

    TEST_CASE("[Faunus] ranunit_polar") {
        Random r;
        int n=2e5;
        Point rtp(0,0,0);
        for (int i=0; i<n; i++)
            rtp += xyz2rtp( ranunit_polar(r) );
        rtp = rtp / n;
        CHECK( rtp.x() == doctest::Approx(1) );
        CHECK( rtp.y() == doctest::Approx(0).epsilon(0.005) ); // theta [-pi:pi] --> <theta>=0
        CHECK( rtp.z() == doctest::Approx(pc::pi/2).epsilon(0.005) );// phi [0:pi] --> <phi>=pi/2
    }
#endif

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
            virtual bool collision(const Point&) const=0; //!< Overlap with boundaries
            virtual Point getLength() const=0; //!< Side lengths

            inline double sqdist( const Point &a, const Point &b ) const {
                return vdist(a,b).squaredNorm();
            } //!< Squared (minimum) distance between two points

            std::string name;

            virtual ~GeometryBase();

            inline BoundaryFunction getBoundaryFunc() const {
                return [this](Point &i){boundary(i);};
            } //!< returns lambda to boundary()

            inline DistanceFunction getDistanceFunc() const {
                return [this](const Point &i, const Point &j){return vdist(i,j);};
            } //!< returns lambda to vdist()

        }; //!< Base class for all geometries

        /**
         * @brief Geometry class for spheres, cylinders, cuboids, hexagonal prism, truncated octahedron, slits, 2D hypersphere
         *
         * All geometries shares the same distance calculation where
         * PBC is disabled by artificially setting long sidelengths.
         *
         * @note
         * - [Efficient Coding of the Minimum Image Convention](http://doi.org/kvs)
         * - [Fast Coding of the Minimum Image Convention](http://doi.org/ck2nrd)
         *
         * @todo Implement unit tests
         */
        class Chameleon : public GeometryBase {
            public:
                enum Variant {CUBOID=0, SPHERE, CYLINDER, SLIT, HEXAGONAL, OCTAHEDRON, HYPERSPHERE2D};
                Variant type;
            private:
                /**
                 * `pbc_disable` is used to disable PBC in `sphere`, `cylinder`, `hexagonal`, `octahedron`, `slit`, `hypersphere2d` etc.
                 * by scaling the (internal) box length by a large number.
                 */
                static constexpr double pbc_disable=100;
                static const std::map<std::string,Variant> names; //!< geometry names
                double radius=0, c1, c2;
                Point len, len_half, len_inv;

                template<typename T=double>
                    inline int anint( T x ) const {
                        return int(x > 0.0 ? x + 0.5 : x - 0.5);
                    } //!< Round to int

            public:
                void setLength(const Point &l);
                double getVolume(int dim=3) const override;
                double getRadius() const;
                Point setVolume(double V, VolumeMethod method=ISOTROPIC) override;
                Point getLength() const override; //!< Enscribed box
                void randompos( Point &m, Random &rand ) const override;
                bool collision(const Point &a) const override;
                void from_json(const json &j);
                void to_json(json &j) const;

                /**
                 * @todo 'HYPERSPHERE2D' does not have a boundary. Does that matter for the program?
                 */
                inline void boundary( Point &a ) const override {
                    if ( type!=HEXAGONAL and type!=OCTAHEDRON and type!=HYPERSPHERE2D ) {
                        if ( std::fabs(a.x()) > len_half.x())
                            a.x() -= len.x() * anint(a.x() * len_inv.x());
                        if ( std::fabs(a.y()) > len_half.y())
                            a.y() -= len.y() * anint(a.y() * len_inv.y());
                        if ( std::fabs(a.z()) > len_half.z())
                            a.z() -= len.z() * anint(a.z() * len_inv.z());
                    } else if ( type == HEXAGONAL ) {
                        const double sqrtThreeByTwo = 0.5*sqrt(3.0);
                        const Point unitvX = {1,0,0};
                        const Point unitvY = {0.5,sqrtThreeByTwo,0};
                        const Point unitvZ = {-0.5,sqrtThreeByTwo,0.0};

                        double tmp = a.dot(unitvX);
                        if (std::fabs(tmp) > len_half.x())
                            a -= len.x() * anint(tmp * len_inv.x())*unitvX;

                        if (a.dot(unitvY) > len_half.x()) {
                            a = a - len.x()*unitvY;
                            if (a.dot(unitvX) < -len_half.x()) // Check that point did not get past x-limit
                                a = a + len.x()*unitvX;
                        }
                        if (a.dot(unitvY) < -len_half.x()) {
                            a = a + len.x()*unitvY;
                            if (a.dot(unitvX) > len_half.x()) // Check that point did not get past x-limit
                                a = a - len.x()*unitvX;
                        }

                        tmp = a.dot(unitvZ);
                        if (std::fabs(tmp) > len_half.x())
                            a -= len.x() * anint(tmp * len_inv.x())*unitvZ;
                        if (std::fabs(a.z()) > len_half.z())
                            a.z() -= len.z() * anint(a.z() * len_inv.z());

                    } else if ( type == OCTAHEDRON ) {
                        const double sqrtThreeI = 1.0/std::sqrt(3.0);
                        const Point unitvXYZ  = Point(1, 1, 1)*sqrtThreeI;
                        const Point unitvXiYZ = Point(1, 1,-1)*sqrtThreeI;
                        const Point unitvXYiZ = Point(1,-1,-1)*sqrtThreeI;
                        const Point unitvXYZi = Point(1,-1, 1)*sqrtThreeI;

                        bool outside = false;
                        do {
                            outside = false;
                            double tmp = a.dot(unitvXYZ);
                            if ( std::fabs(tmp) > len_half.x()) {
                                a -= len.x() * anint(tmp * len_inv.x()) * unitvXYZ;
                                outside = true;
                            }
                            tmp = a.dot(unitvXiYZ);
                            if ( std::fabs(tmp) > len_half.x()) {
                                a -= len.x() * anint(tmp * len_inv.x()) * unitvXiYZ;
                                outside = true;
                            }
                            tmp = a.dot(unitvXYiZ);
                            if ( std::fabs(tmp) > len_half.x()) {
                                a -= len.x() * anint(tmp * len_inv.x()) * unitvXYiZ;
                                outside = true;
                            }
                            tmp = a.dot(unitvXYZi);
                            if ( std::fabs(tmp) > len_half.x()) {
                                a -= len.x() * anint(tmp * len_inv.x()) * unitvXYZi;
                                outside = true;
                            }
                        } while (outside);

                        if ( std::fabs(a.x()) > len_half.y())
                            a.x() -= len.y() * anint(a.x() * len_inv.y());

                        if ( std::fabs(a.y()) > len_half.y())
                            a.y() -= len.y() * anint(a.y() * len_inv.y());

                        if ( std::fabs(a.z()) > len_half.y())
                            a.z() -= len.y() * anint(a.z() * len_inv.y());
                    }
                } //!< Apply boundary conditions

                inline Point vdist(const Point &a, const Point &b) const override {
                    Point r(a - b);
                    if (type!=HEXAGONAL and type!=OCTAHEDRON and type!=HYPERSPHERE2D) {
                        if ( r.x() > len_half.x())
                            r.x() -= len.x();
                        else if ( r.x() < -len_half.x())
                            r.x() += len.x();

                        if ( r.y() > len_half.y())
                            r.y() -= len.y();
                        else if ( r.y() < -len_half.y())
                            r.y() += len.y();

                        if ( r.z() > len_half.z())
                            r.z() -= len.z();
                        else if ( r.z() < -len_half.z())
                            r.z() += len.z();
                    } else if(type == HYPERSPHERE2D) { // ugly but works, needs fixing though...
		      double angle = std::acos(a.dot(b)/radius/radius);
		      double distance = radius*angle;
		      return r/r.norm()*distance;
		    } else
                        boundary(r);
                    return r;
                }
        };

        void to_json(json&, const Chameleon&);
        void from_json(const json&, Chameleon&);

#ifdef DOCTEST_LIBRARY_INCLUDED
        TEST_CASE("[Faunus] Chameleon") {

            using doctest::Approx;
            Chameleon geo;
            Random slump;

            SUBCASE("cuboid") {
                geo = R"( {"type": "cuboid", "length": [2,3,4]} )"_json;
                CHECK( geo.getVolume() == doctest::Approx(2*3*4) );
                Point a(1.1, 1.5, -2.001);
                CHECK( geo.collision(a) == true);
                geo.getBoundaryFunc()(a);
                CHECK( geo.collision(a) == false);
                CHECK( a.x() == Approx(-0.9) );
                CHECK( a.y() == Approx(1.5) );
                CHECK( a.z() == Approx(1.999) );
                CHECK( geo.vdist({1,2,3}, a) == geo.getDistanceFunc()({1,2,3},a) );
                Point b = a;
                geo.boundary(b);
                CHECK( a == b );

                // check that geometry is properly enscribed in a cuboid
                Point L = geo.getLength();
                CHECK( L.x() == Approx(2) );
                CHECK( L.y() == Approx(3) );
                CHECK( L.z() == Approx(4) );

                // check random position
                bool containerOverlap=false;
                for (int i=0; i<1e4; i++) {
                    geo.randompos(a, slump);
                    if (geo.collision(a))
                        containerOverlap=true;
                }
                CHECK( containerOverlap==false );
            }

            SUBCASE("slit") {
                geo = R"( {"type": "slit", "length": [2,3,4]} )"_json;
                CHECK( geo.getVolume() == doctest::Approx(2*3*4) );
                Point a(1.1, 1.5, -2.001);
                CHECK( geo.collision(a) == true);
                geo.getBoundaryFunc()(a);
                CHECK( geo.collision(a) == true);
                CHECK( a.x() == doctest::Approx(-0.9) );
                CHECK( a.y() == doctest::Approx(1.5) );
                CHECK( a.z() == doctest::Approx(-2.001) );

                Point b = a;
                b.z() = -b.z();
                geo.boundary(b);
                CHECK( b.z() == doctest::Approx(2.001) );
                CHECK( geo.vdist(a,b).x()==0);
                CHECK( geo.vdist(a,b).y()==0);
                CHECK( geo.vdist(a,b).z()==-4.002);

                // check that geometry is properly enscribed in a cuboid
                Point L = geo.getLength();
                CHECK( L.x() == Approx(2) );
                CHECK( L.y() == Approx(3) );
                CHECK( L.z() == Approx(4) );

                // check random position
                bool containerOverlap=false;
                for (int i=0; i<1e4; i++) {
                    geo.randompos(a, slump);
                    if (geo.collision(a))
                        containerOverlap=true;
                }
                CHECK( containerOverlap==false );
            }

            SUBCASE("cylinder") {
                json j = {  {"type","cylinder"}, {"radius",1.0}, {"length", 1/pc::pi} };
                geo = j;
                CHECK( geo.getVolume() == doctest::Approx( 1.0 ) );
                CHECK( geo.collision({1.01,0,0}) == true);
                CHECK( geo.collision({0.99,0,0}) == false);
                CHECK( geo.collision({0.99,1/pc::pi+0.01,0}) == true);

                // check that geometry is properly enscribed in a cuboid
                Point L = geo.getLength();
                CHECK( L.x() == Approx(2) );
                CHECK( L.y() == Approx(2) );
                CHECK( L.z() == Approx(1/pc::pi) );

                // check random position
                Point a;
                bool containerOverlap=false;
                for (int i=0; i<1e4; i++) {
                    geo.randompos(a, slump);
                    if (geo.collision(a))
                        containerOverlap=true;
                }
                CHECK( containerOverlap==false );
            }

            SUBCASE("hexagonal") {
                json j = {  {"type","hexagonal"}, {"radius",1.0}, {"length", 1.0/2.0/std::sqrt(3.0)} };
                geo = j;
                CHECK( geo.getVolume() == doctest::Approx( 1.0 ) );
                CHECK( geo.collision({1.01,0,0}) == true);
                CHECK( geo.collision({0.99,0,0}) == false);
                CHECK( geo.collision({0.0,2.0/std::sqrt(3.0)+0.01,0}) == true);

                // check that geometry is properly enscribed in a cuboid
                Point L = geo.getLength();
                CHECK( L.x() == Approx(1.0) );
                CHECK( L.y() == Approx(1.0) );
                CHECK( L.z() == Approx(1.0/2.0/std::sqrt(3.0)) );

                // check random position
                Point a;
                bool containerOverlap=false;
                for (int i=0; i<1e4; i++) {
                    geo.randompos(a, slump);
                    if (geo.collision(a))
                        containerOverlap=true;
                }
                CHECK( containerOverlap==false );
            }

            SUBCASE("sphere") {
                geo = R"( { "type": "sphere", "radius": 2 } )"_json;
                CHECK( geo.getVolume() == doctest::Approx( 4*pc::pi*8/3.0 ) );
                CHECK( geo.collision( { 2.01, 0, 0 } ) == true );
                CHECK( geo.collision( { 1.99, 0, 0 } ) == false );
                CHECK( geo.getLength().squaredNorm() == doctest::Approx(48) );

                // check that geometry is properly enscribed in a cuboid
                Point L = geo.getLength();
                CHECK( L.x() == Approx(4) );
                CHECK( L.y() == Approx(4) );
                CHECK( L.z() == Approx(4) );

                geo.setVolume(123.4);
                CHECK( geo.getVolume() == doctest::Approx(123.4) );

                // check random position
                Point a;
                bool containerOverlap=false;
                for (int i=0; i<1e4; i++) {
                    geo.randompos(a, slump);
                    if (geo.collision(a))
                        containerOverlap=true;
                }
                CHECK( containerOverlap==false );
            }
        }
#endif

        /*
           void unwrap( Point &a, const Point &ref ) const {
           a = vdist(a, ref) + ref;
           } //!< Remove PBC with respect to a reference point
         */
        enum class weight { MASS, CHARGE, GEOMETRIC };

        template<typename Titer, typename Tparticle=typename Titer::value_type, typename weightFunc>
            Point anyCenter( Titer begin, Titer end, BoundaryFunction boundary, const weightFunc &weight,
                    const Point &shift={0,0,0})
            {
                double sum=0;
                Point c(0,0,0), t;
                for (auto &i=begin; i!=end; ++i) {
                    t = i->pos + shift;       // translate to origo
                    boundary(t);
                    double w = weight(*i);
                    c += w * t;
                    sum += w;
                }
                if (sum<=pc::epsilon_dbl)
                    return {0,0,0};
                else
                    c = c/sum - shift;
                boundary(c);
                return c;
            } //!< Mass, charge, or geometric center of a collection of particles

        template<typename Titer, typename Tparticle=typename Titer::value_type>
            Point massCenter(Titer begin, Titer end, BoundaryFunction boundary=[](Point&){}, const Point &shift={0,0,0}) {
                return anyCenter(begin, end, boundary,
                        []( const Tparticle &p ){ return atoms.at(p.id).mw; }, shift );
            } // Mass center

#ifdef DOCTEST_LIBRARY_INCLUDED
        TEST_CASE("[Faunus] anyCenter") {
            typedef Particle<Radius, Charge, Dipole, Cigar> T;

            Chameleon cyl = json( {{"type","cuboid"}, {"length",100}, {"radius",20}} );
            std::vector<T> p;

            CHECK( !atoms.empty() ); // set in a previous test
            p.push_back( atoms[0] );
            p.push_back( atoms[0] );

            p.front().pos = {10, 10, -10};
            p.back().pos  = {15, -10, 10};

            Point cm = Geometry::massCenter(p.begin(), p.end(), cyl.getBoundaryFunc());
            CHECK( cm.x() == doctest::Approx(12.5) );
            CHECK( cm.y() == doctest::Approx(0) );
            CHECK( cm.z() == doctest::Approx(0) );
        }
#endif

        template<class Titer=typename std::vector<T>::iterator>
            void translate( Titer begin, Titer end, const Point &d,
                    BoundaryFunction boundary=[](Point&){}  )
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
