#pragma once

#include "core.h"

/** @brief Faunus main namespace */
namespace Faunus {

    /** @brief Simulation geometries and related operations */
    namespace Geometry {

        typedef std::function<void(Point&)> BoundaryFunction;
        typedef std::function<void(const Point&, const Point&)> DistanceFunction;

        struct GeometryBase {
            virtual void setVolume(double, const std::vector<double>&)=0; //!< Set volume
            virtual double getVolume(int=3) const=0; //!< Get volume
            virtual void randompos( Point&, Random& ) const=0; //!< Generate random position
            virtual Point vdist( const Point&, const Point& ) const=0; //!< (Minimum) distance between two points
            virtual void boundary( Point& ) const=0; //!< Apply boundary conditions
            virtual bool collision( const Point&, double=0) const=0; //!< Check if position lies within

            BoundaryFunction boundaryFunc; //!< Functor for boundary()
            DistanceFunction distanceFunc; //!< Functor for vdist()

            std::string name;

            GeometryBase() {
                using namespace std::placeholders;
                boundaryFunc = std::bind( &GeometryBase::boundary, this, _1);
                distanceFunc = std::bind( &GeometryBase::vdist, this, _1, _2);
            }

        }; //!< Base class for all geometries

        /** @brief Cuboidal box */
        class Box : public GeometryBase {
            protected:
                Point len, //!< side length
                      len_half, //!< half side length
                      len_inv; //!< inverse side length

            public:

                void setLength( const Point &l ) {
                    len = l;
                    len_half = l*0.5;
                    len_inv = l.cwiseInverse();
                } //!< Set cuboid side length

                void setVolume(double V, const std::vector<double> &s={1,1,1}) override {
                    double l = std::cbrt(V);
                    setLength( {l,l,l} );
                }

                double getVolume(int dim=3) const override {
                    assert(dim==3);
                    return len.x() * len.y() * len.z();
                }

                const Point& getLength() const { return len; } //!< Side lengths

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
                    if (b.getVolume()<=0)
                        throw std::runtime_error("volume is zero or less");
                }
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

                void setVolume(double V, const std::vector<double> &s={}) override {
                    r = std::sqrt( V / (pc::pi * len) );
                    set( r, len );
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

                bool collision(const Point &a, double r) const override
                {
                    return (a.squaredNorm() > r2) ? true : false;
                }

        };

        enum class weight { MASS, CHARGE, GEOMETRIC };

        template<typename Titer, typename Tparticle=typename Titer::value_type, typename weightFunc, typename boundaryFunc>
            Point anyCenter( Titer begin, Titer end, const boundaryFunc &boundary, const weightFunc &weight)
            {
                double sum=0;
                Point c(0,0,0);
                try {
                    for (auto &i=begin; i!=end; ++i) {
                        double w = weight(*i);
                        c += w * i->pos;
                        sum += w;
                    }
                    c /= sum;
                    if (std::isnan(c[0]))
                        throw std::runtime_error("sum of weights is zero");
                }
                catch(std::exception& e) {
                    throw std::runtime_error("anyCenter: " + std::string(e.what()));
                }
                return c;
            } //!< Mass, charge, or geometric center of a collection of particles

        template<typename Titer, typename Tparticle=typename Titer::value_type>
            Point massCenter(Titer begin, Titer end, const BoundaryFunction &boundary=[](Point&){}) {
                return anyCenter(begin, end, boundary,
                        []( const Tparticle &p ){ return atoms<Tparticle>.at(p.id).mw; } );
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
                    const BoundaryFunction &boundary=[](Point &i){}  )
            {
                for ( auto i=begin; i!=end; ++i )
                {
                    i->pos += d;
                    boundary(i->pos);
                }
            } //!< Vector displacement of a range of particles

        template<typename Titer>
            void cm2origo( Titer begin, Titer end, const BoundaryFunction &boundary=[](Point&){} )
            {
                Point cm = massCenter(begin, end, boundary);
                translate(begin, end, -cm, boundary);
            } //!< Translate to that mass center is in (0,0,0)


        template<typename Titer>
            void rotate(
                    Titer begin,
                    Titer end,
                    const Eigen::Quaterniond &q,
                    const BoundaryFunction &boundary=[](Point&){} ,
                    const Point &shift=Point(0,0,0) )
            {
                auto m = q.toRotationMatrix(); // rotation matrix
                for (auto &i=begin; i!=end; ++i) {
                    i->rotate(q,m); // rotate internal coordinates
                    i->pos += shift;
                    boundary(i->pos);
                    i->pos = q*i->pos - shift;
                    boundary(i->pos);
                }
            } //!< Rotate particle pos and internal coordinates

    } //geometry namespace
}//end of faunus namespace
