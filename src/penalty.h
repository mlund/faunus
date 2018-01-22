#pragma once

#include "core.h"
#include "group.h"
#include "space.h"

namespace Faunus {

    namespace ReactionCoordinate {

        /**
         * @brief Base class for reaction coordinates
         *
         * Derived classes must implement the function operator.
         *
         * Keyword  | Description
         * :------- | :----------------
         * `dim`    | Number of dimensions
         * `min`    | minumum coordinate values (array)
         * `max`    | maximum coordinate values (array)
         * `scale`  | penalty scaling (array)
         * `nstep`  | update frequency for penalty updates, sampling etc. (steps)
         */
        class ReactionCoordinateBase {
            public:
                typedef std::vector<double> Tvec;
            protected:
                std::function<Tvec()> f=nullptr; // returns reaction coordinate
            public:
                inline virtual void _to_json(json &j) const {};
                size_t nstep;     // update frequency [steps]
                Tvec min, max, scale;
                std::string name;

                size_t dim() const { return min.size(); }   //!< dimension of coordinate

                Tvec operator()() {
                    assert(f!=nullptr);
                    return f();
                }

                bool inrange( const Tvec &coord ) const
                {
                    for ( size_t i = 0; i < coord.size(); ++i )
                        if ( coord[i] < min[i] || coord[i] > max[i] )
                            return false;
                    return true;
                } //!< Determines if coordinate is within [min,max]
        };

        void to_json(json &j, const ReactionCoordinateBase &r) {
            j = {{"dim", r.dim()}, {"min", r.min}, {"max", r.max}};
            if (r.nstep>0) {
                j["nstep"] = r.nstep;
                j["scale"] = r.scale;
            }
            r._to_json(j);
        }

        void from_json(const json &j, ReactionCoordinateBase &r) {
            size_t d = j.value("dim", 1);
            r.min.resize(d, -pc::infty);
            r.min = j.value("min", r.min);
            r.max.resize(d, pc::infty);
            r.max = j.value("max", r.max);
            r.scale.resize(d, 1);
            r.scale = j.value("scale", r.scale);
            r.nstep = j.value("nstep", 0);

            if ( (r.min.size()!=r.max.size() ) || r.min.size()!=r.scale.size()
                    || r.min.size()<1 || d!=r.dim())
                throw std::runtime_error(
                        "Reaction coordinate error: min, max, scale must have equal length >=1");
        }

#ifdef DOCTEST_LIBRARY_INCLUDED
        TEST_CASE("[Faunus] ReactionCoordinateBase")
        {
            using doctest::Approx;
            ReactionCoordinateBase c = R"({"dim":2, "min":[0,-1.5], "max":[8.5,7], "nstep":10})"_json;
            CHECK( c.dim() == 2 );
            CHECK( c.min.size() == 2 );
            CHECK( c.max.size() == 2 );
            CHECK( c.scale.size() == 2 );
            CHECK( c.min[0] == 0 );
            CHECK( c.min[1] == -1.5 );
            CHECK( c.max[0] == 8.5 );
            CHECK( c.max[1] == 7 );
            CHECK( c.nstep == 10) ;
        }
#endif

        /**
         * @brief Reaction coordinate: molecule-molecule mass-center separation
         */
        class MassCenterSeparation : public ReactionCoordinateBase {
            public:
                Point dir={1,1,1};
                int id1, id2;
                void _to_json(json &j) const override {
                    j["dir"] = dir;
                }
                template<class Tspace>
                    MassCenterSeparation( const json &j, Tspace &spc ) {
                        from_json(j, *this);
                        dir = j.value("dir", dir);
                        auto i = j.value("index", std::vector<int>({0,0}) );
                        if (i.size()!=2 || dir.size()!=3)
                            throw std::runtime_error("CM reaction coordinate error");

                        //string name1 = j.at("first");
                        //string name2 = j.at("second");
                        //g1 = s.findMolecules(name1, true).at(i[0]);
                        //g2 = s.findMolecules(name2, true).at(i[1]);

                        //assert(!g1->empty());
                        //assert(!g2->empty());
                        f = [&spc, dir=dir]() {
                            Point cm1;// = massCenter(s.geo, p, *g1);
                            Point cm2;// = massCenter(s.geo, p, *g2);
                            return Tvec({spc.geo.vdist(cm1, cm2).cwiseProduct(dir).norm()});
                        };
                    }
        };
#ifdef DOCTEST_LIBRARY_INCLUDED
        TEST_CASE("[Faunus] MassCenterSeparation")
        {
            using doctest::Approx;
            typedef Space<Geometry::Cuboid, Particle<>> Tspace;
            Tspace spc;
            MassCenterSeparation c( R"({"dir":[1,1,0]})"_json, spc);
            CHECK( c.dim() == 1 );
            CHECK( c.min[0] == -pc::infty );
            CHECK( c.max[0] == pc::infty );
            CHECK( c.nstep == 0 );
            CHECK( c.dir.x() == 1 );
            CHECK( c.dir.y() == 1 );
            CHECK( c.dir.z() == 0 );
        }
#endif

    } // namespace
} // namespace
