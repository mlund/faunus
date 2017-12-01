#pragma once

#include "core.h"
#include "auxiliary.h"

namespace Faunus {
    namespace Potential {

        struct PairPotentialBase {
            std::string name;
            PairPotentialBase() {}
            PairPotentialBase(const std::string &name) : name(name) {}
            virtual void from_json(const json&)=0;
        }; //!< Base for pair-potentials

        void from_json(const json &j, PairPotentialBase &pairpot) {
            pairpot.from_json(j);
        } //!< Configure via json

        template<class T1, class T2>
            struct CombinedPairPotential : public PairPotentialBase {
                T1 first;  //!< First pair potential of type T1
                T2 second; //!< Second pair potential of type T2

                CombinedPairPotential(const T1 &a, const T2 &b) : first(a), second(b) {}

                template<typename... T>
                    double operator()(const Particle<T...> &a, const Particle<T...> &b, const Point &r) const {
                        return first(a, b, r) + second(a, b, r);
                    }

                void from_json(const json &j) override {
                    first = j;
                    second = j;
                }
            };

        template<class T1, class T2,
            class = typename std::enable_if<std::is_base_of<PairPotentialBase, T1>::value>::type,
            class = typename std::enable_if<std::is_base_of<PairPotentialBase, T2>::value>::type>
                CombinedPairPotential<T1, T2> &operator+(const T1 &pot1, const T2 &pot2) {
                    return *(new CombinedPairPotential<T1, T2>(pot1, pot2));
                } //!< Add two pair potentials

        struct Dummy : public PairPotentialBase {
            template<typename... T>
                double operator()(const Particle<T...> &a, const Particle<T...> &b, const Point &r) const {
                    return 0;
                }
            void from_json(const json &) override {}
        }; //!< A dummy pair potential that always returns zero

        template<typename Tparticle>
            struct SigmaEpsilonTable {
                enum Mixers {LB};
                Mixers mixer = LB;
                PairMatrix<double> s2,eps; // matrix of sigma_ij^2 and 4*eps_ij
                //decltype(atoms<Tparticle>)&  = atoms<Tparticle>;
            }; //!< Table of sigma and epsilons

        template<typename Tparticle>
            void from_json(const json &j, SigmaEpsilonTable<Tparticle> &m) {
                std::function<std::pair<double,double>(double,double,double,double)> mixerFunc;

                std::string mixer = j.value("mixer", std::string("LB"));
                if (mixer=="LB")
                    m.mixer=SigmaEpsilonTable<Tparticle>::LB;
                switch(m.mixer) {
                    case SigmaEpsilonTable<Tparticle>::LB:
                        mixerFunc = [](double s1, double s2, double e1, double e2) {
                            return std::pair<double,double>( { (s1+s2)/2, std::sqrt(e1*e2) } );
                        };
                        break;
                    default:
                        throw std::runtime_error("unknown mixing rule");
                }
                size_t n=atoms<Tparticle>.size(); // number of atom types
                m.s2.resize(n); // not required...
                m.eps.resize(n);// ...but possible reduced mem. fragmentation
                for (auto &i : atoms<Tparticle>)
                    for (auto &j : atoms<Tparticle>) {
                        double sigma, epsilon; // mixed values
                        std::tie( sigma, epsilon ) = mixerFunc(i.sigma, j.sigma, i.eps, j.eps);
                        m.s2.set(  i.id(), j.id(), sigma*sigma );
                        m.eps.set( i.id(), j.id(), 4*epsilon ); // should already be in kT
                    }

                // custom eps/sigma for specific pairs
                auto &_j = j["ljcustom"];
                if (_j.is_object()) {
                    for (auto it=_j.begin(); it!=_j.end(); ++it) {
                        auto v = words2vec<std::string>( it.key() );
                        if (v.size()==2) {
                            int id1 = (*findName( atoms<Tparticle>, v[0])).id();
                            int id2 = (*findName( atoms<Tparticle>, v[1])).id();
                            m.s2.set( id1, id2, std::pow( it.value().at("sigma").get<double>(), 2) );
                            m.eps.set(id1, id2, 4*it.value().at("eps").get<double>() * 1.0_kJmol);
                        } else
                            std::runtime_error("custom LJ parameters require exactly two atoms");
                    }
                }
            }

        /**
         * @brief Lennard-Jones with arbitrary mixing rule
         */
        template<typename Tparticle>
            struct LennardJones : public PairPotentialBase {
                LennardJones() { name="lennardjones"; }
                SigmaEpsilonTable<Tparticle> m; // table w. sigma_ij^2 and 4xepsilon

                template<typename... T>
                    Point force(const Particle<T...> &a, const Particle<T...> &b, double r2, const Point &p) const {
                        double s6=powi<3>( m.s2(a.id,b.id) );
                        double r6=r2*r2*r2;
                        double r14=r6*r6*r2;
                        return 6.*m.eps(a.id,b.id) * s6 * (2*s6-r6) / r14 * p;
                    }

                template<typename... T>
                    double operator()(const Particle<T...> &a, const Particle<T...> &b, const Point &r) const {
                        double x=m.s2(a.id,b.id)/r.squaredNorm(); //s2/r2
                        x=x*x*x; // s6/r6
                        return m.eps(a.id,b.id) * (x*x - x);
                    }

                void from_json(const json &j) override { m = j.at(name); }
            };

        struct Coulomb : public PairPotentialBase {
            Coulomb() { name="coulomb"; }
            double lB; //!< Bjerrum length
            template<typename... T>
                double operator()(const Particle<T...> &a, const Particle<T...> &b, const Point &r) const {
                    return lB * a.charge * b.charge / r.norm();
                }
            void from_json(const json &j) override {
                lB = pc::lB( j.at(name).at("epsr").get<double>() );
            }
            void to_json(json &j) {
                j[name] = {"lB", lB};
            }
        };

        struct HardSphere : public PairPotentialBase {
            template<typename... T>
                double operator()(const Particle<T...> &a, const Particle<T...> &b, const Point &r2) const {
                    auto d=a.radius+b.radius;
                    return (r2.squaredNorm()<d*d) ? pc::infty : 0;
                }
            void from_json(const json &j) override {}
        };
    }//end of namespace Potential
}//end of namespace Faunus
