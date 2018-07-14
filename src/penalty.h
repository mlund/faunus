#pragma once

#include "core.h"
#include "group.h"
#include "space.h"
#include <Eigen/Dense>

namespace Faunus {

    namespace ReactionCoordinate {

        /**
         * @brief Base class for reaction coordinates
         */
        struct ReactionCoordinateBase {
            std::function<double()> f=nullptr; // returns reaction coordinate
            inline virtual void _to_json(json &j) const {};
            inline virtual double normalize(double) const { return 1.; }
            double binwidth=0, min=0, max=0;
            std::string name;

            double operator()() {
                assert(f!=nullptr);
                return f();
            } //!< Calculates reaction coordinate

            bool inRange(double coord) const {
                return (coord>=min && coord<=max);
            } //!< Determines if coordinate is within [min,max]
        };

        void to_json(json &j, const ReactionCoordinateBase &r) {
            j = { {"range",{r.min,r.max}}, {"resolution",r.binwidth} };
            r._to_json(j);
        }

        void from_json(const json &j, ReactionCoordinateBase &r) {
            r.binwidth = j.value("resolution", 0.5);
            auto range = j.value("range", std::vector<double>({0,0}));
            if (range.size()==2)
                if (range[0]<=range[1]) {
                    r.min = range[0];
                    r.max = range[1];
                    return;
                 }
            throw std::runtime_error(r.name + ": 'range' require two numbers: [min, max>=min]");
        }

#ifdef DOCTEST_LIBRARY_INCLUDED
        TEST_CASE("[Faunus] ReactionCoordinateBase")
        {
            using doctest::Approx;
            ReactionCoordinateBase c = R"({"range":[-1.5, 2.1], "resolution":0.2})"_json;
            CHECK( c.min == Approx(-1.5) );
            CHECK( c.max == Approx( 2.1) );
            CHECK( c.binwidth == Approx(0.2) );
            CHECK( c.inRange(-1.5)  == true );
            CHECK( c.inRange(-1.51) == false );
            CHECK( c.inRange( 2.11) == false );
            CHECK( c.inRange( 2.1)  == true );
        }
#endif

        class SystemProperty : public ReactionCoordinateBase {
            protected:
                std::string property;
            public:
                template<class Tspace>
                    SystemProperty(const json &j, Tspace &spc) {
                        name = "system";
                        from_json(j, *this);
                        property = j.at("property").get<std::string>();
                        if (property=="V") f = [&g=spc.geo]() { return g.getVolume();};
                        if (property=="height") f = [&g=spc.geo]() { return g.getLength().z(); };
                        if (f==nullptr)
                            throw std::runtime_error(name + ": unknown property '" + property + "'");
                    }
                void _to_json(json &j) const override {
                    j["property"] = property;
                }
        };

        class AtomProperty : public ReactionCoordinateBase {
            protected:
                std::string property;
                size_t index; // atom index
            public:
                template<class Tspace>
                    AtomProperty(const json &j, Tspace &spc) {
                        name = "atom";
                        from_json(j, *this);
                        index = j.at("index");
                        property = j.at("property").get<std::string>();
                        if (property=="x") f = [&i=spc.p.at(index).pos]() { return i.x(); };
                        if (property=="y") f = [&i=spc.p.at(index).pos]() { return i.y(); };
                        if (property=="z") f = [&i=spc.p.at(index).pos]() { return i.z(); };
                        if (property=="R") f = [&i=spc.p.at(index).pos]() { return i.norm(); };
                        if (f==nullptr)
                            throw std::runtime_error(name + ": unknown property '" + property + "'");
                    }
                void _to_json(json &j) const override {
                    j["property"] = property;
                    j["index"] = index;
                }
        };

        /**
         * @brief Reaction coordinate: molecule-molecule mass-center separation
         */
        struct MassCenterSeparation : public ReactionCoordinateBase {
            Eigen::Vector3i dir={1,1,1};
            std::vector<size_t> index;
            std::vector<std::string> type;
            template<class Tspace>
                MassCenterSeparation(const json &j, Tspace &spc) {
                    typedef typename Tspace::Tparticle Tparticle;
                    name = "cmcm";
                    from_json(j, *this);
                    dir = j.value("dir", dir);
                    index = j.at("index").get<decltype(index)>();
                    type = j.at("type").get<decltype(type)>();
                    if (index.size()==2) {
                        f = [&spc, dir=dir, i=index[0], j=index[1]]() {
                            auto &cm1 = spc.groups[i].cm;
                            auto &cm2 = spc.groups[j].cm;
                            return spc.geo.vdist(cm1, cm2).cwiseProduct(dir.cast<double>()).norm(); 
                        };
                    }
                    else if (type.size()==2) {
                        f = [&spc, dir=dir, type1=type[0], type2=type[1]]() {
                            Group<Tparticle> g(spc.p.begin(), spc.p.end());
                            auto slice1 = g.find_id(findName(atoms<Tparticle>, type1)->id());
                            auto slice2 = g.find_id(findName(atoms<Tparticle>, type2)->id());
                            auto cm1 = Geometry::massCenter(slice1.begin(), slice1.end(), spc.geo.boundaryFunc);
                            auto cm2 = Geometry::massCenter(slice2.begin(), slice2.end(), spc.geo.boundaryFunc);
                            return spc.geo.vdist(cm1, cm2).cwiseProduct(dir.cast<double>()).sum();
                        };
                    }
                    else
                        throw std::runtime_error(name + ": specify exactly two molecule indeces or two atom types");
                }

            double normalize(double coord) const override {
                int dim=dir.sum();
                if (dim==2) return 1/(2*pc::pi*coord);
                if (dim==3) return 1/(4*pc::pi*coord*coord);
                return 1.0;
            } // normalize by volume element

            void _to_json(json &j) const override {
                j["dir"] = dir;
                j["index"] = index;
                j["type"] = type;
            }
        };
        /**
         * @brief Reaction coordinate: angle between principal axis and cartesian axis
         */
        struct PrincipalAxisAngle : public ReactionCoordinateBase {
            Eigen::Vector3d dir={0,0,1};
            size_t index;
            template<class Tspace>
                PrincipalAxisAngle(const json &j, Tspace &spc) {
                    typedef typename Tspace::Tparticle Tparticle;
                    name = "angle";
                    from_json(j, *this);
                    dir = j.value("dir", dir);
                    index = j.at("index").get<decltype(index)>();
                    auto name = molecules<decltype(spc.p)>.at(spc.groups[index].id).name;
                    cout << "Molecule Name: " << name << endl;
                    f = [&spc, dir=dir, i=index]() {
                        auto &cm = spc.groups[i].cm;
                        //Point vec = spc.geo.vdist(spc.groups[i].begin()->pos,(spc.groups[i].end()-1)->pos);
                        //vec = vec / vec.norm();
                        //cout << "P1 " << atoms<Tparticle>[spc.groups[i].begin()->id].name << endl;
                        //cout << "P2 " << atoms<Tparticle>[(spc.groups[i].end()-1)->id].name << endl;
                        auto S = Geometry::gyration(spc.groups[i].begin(), spc.groups[i].end(), spc.geo.boundaryFunc, cm);
                        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> esf(S);
                        Point eivals = esf.eigenvalues();
                        std::ptrdiff_t i_eival;
                        double minOfeivals = eivals.minCoeff(&i_eival);
                        Point vec = esf.eigenvectors().col(i_eival).real();
                        double cosine = vec.dot(dir);
                        double angle = acos(abs(cosine)) * 180. / pc::pi;
                        return angle; 
                    };
                }

            void _to_json(json &j) const override {
                j["dir"] = dir;
                j["index"] = index;
            }
        };

#ifdef DOCTEST_LIBRARY_INCLUDED
        TEST_CASE("[Faunus] MassCenterSeparation")
        {
            using doctest::Approx;
            typedef Space<Geometry::Cuboid, Particle<>> Tspace;
            Tspace spc;
            MassCenterSeparation c( R"({"dir":[1,1,0], "index":[7,8], "type":[] })"_json, spc);
            CHECK( c.dir.x() == 1 );
            CHECK( c.dir.y() == 1 );
            CHECK( c.dir.z() == 0 );
            CHECK( c.index == decltype(c.index)({7,8}) );
        }
#endif

    } // namespace
} // namespace
