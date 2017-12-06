#pragma once

#include "core.h"
#include "auxiliary.h"

namespace Faunus {
    namespace Potential {

        struct PairPotentialBase {
            std::string name;
            virtual void to_json(json&) const=0;
            virtual void from_json(const json&)=0;
        }; //!< Base for all pair-potentials

        void to_json(json &j, const PairPotentialBase &base) {
            base.name.empty() ? base.to_json(j) : base.to_json(j[base.name]);
        } //!< Serialize any pair potential to json

        void from_json(const json &j, PairPotentialBase &base) {
            base.name.empty() ? base.from_json(j) : base.from_json(j.at(base.name));
        } //!< Serialize any pair potential from json

        template<class T1, class T2>
            struct CombinedPairPotential : public PairPotentialBase {
                T1 first;  //!< First pair potential of type T1
                T2 second; //!< Second pair potential of type T2
                template<typename... T>
                    double operator()(const Particle<T...> &a, const Particle<T...> &b, const Point &r) const {
                        return first(a, b, r) + second(a, b, r);
                    }

                void from_json(const json &j) override {
                    first = j;
                    second = j;
                }

                void to_json(json &j) const override { j = {first,second}; }
            };

        template<class T1, class T2,
            class = typename std::enable_if<std::is_base_of<PairPotentialBase, T1>::value>::type,
            class = typename std::enable_if<std::is_base_of<PairPotentialBase, T2>::value>::type>
                CombinedPairPotential<T1, T2> &operator+(const T1 &pot1, const T2 &pot2) {
                    return *(new CombinedPairPotential<T1, T2>(pot1, pot2));
                } //!< Add two pair potentials

        struct Dummy : public PairPotentialBase {
            Dummy() { name="dummy"; }
            template<typename... T>
                double operator()(const Particle<T...> &a, const Particle<T...> &b, const Point &r) const {
                    return 0;
                }
            void from_json(const json&) override {}
            void to_json(json&) const override {}
        }; //!< A dummy pair potential that always returns zero

        template<typename Tparticle>
            struct SigmaEpsilonTable {
                enum Mixers {LB};
                Mixers mixer = LB;
                PairMatrix<double> s2,eps; // matrix of sigma_ij^2 and 4*eps_ij
            }; //!< Table of sigma and epsilons

        template<typename Tparticle>
            void from_json(const json &j, SigmaEpsilonTable<Tparticle> &m) {
                std::function<std::pair<double,double>(double,double,double,double)> mixerFunc;

                std::string mixer = j.value("mixing", std::string("LB"));
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
                if (j.count("ljcustom")==1) {
                    auto &_j = j.at("ljcustom");
                    if (_j.is_object()) {
                        for (auto it=_j.begin(); it!=_j.end(); ++it) {
                            auto v = words2vec<std::string>( it.key() );
                            if (v.size()==2) {
                                int id1 = (*findName( atoms<Tparticle>, v[0])).id();
                                int id2 = (*findName( atoms<Tparticle>, v[1])).id();
                                m.s2.set( id1, id2, std::pow( it.value().at("sigma").get<double>(), 2) );
                                m.eps.set(id1, id2, 4*it.value().at("eps").get<double>() * 1.0_kJmol);
                            } else
                            std::runtime_error("custom LJ parameters require exactly two space-separated atoms");
                        }
                    }
                }
            }

        template<typename Tparticle>
            void to_json(json &j, const SigmaEpsilonTable<Tparticle> &m) {
                j["mixing"] = "LB = Lorentz-Berthelot";
                j["epsilon unit"] = "kJ/mol";
                auto& _j = j["combinations"];
                for (size_t i=0; i<m.eps.size(); i++)
                    for (size_t j=0; j<m.eps.size(); j++)
                        if (i>=j) {
                            auto str = atoms<Tparticle>[i].name+" "+atoms<Tparticle>[j].name;
                            _j[str] = { {"eps", m.eps(i,j)/4.0_kJmol}, {"sigma", std::sqrt(m.s2(i,j))}  };
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

                void to_json(json &j) const override { j = m; }
                void from_json(const json &j) override { m = j; }
            };

        struct Coulomb : public PairPotentialBase {
            Coulomb() { name="coulomb"; }
            double lB; //!< Bjerrum length
            template<typename... T>
                double operator()(const Particle<T...> &a, const Particle<T...> &b, const Point &r) const {
                    return lB * a.charge * b.charge / r.norm();
                }

            void to_json(json &j) const override { j["lB"] = lB; }
            void from_json(const json &j) override { lB = pc::lB( j.at("epsr") ); }
        };

        template<class Tparticle>
            struct HardSphere : public PairPotentialBase {
                PairMatrix<double> d2; // matrix of (r1+r2)^2
                HardSphere() {
                    name="hardsphere";
                    for (auto &i : atoms<Tparticle>)
                        for (auto &j : atoms<Tparticle>)
                            d2.set( i.id(), j.id(), std::pow((i.sigma+j.sigma)/2,2));
                }
                template<typename... T>
                    double operator()(const Particle<T...> &a, const Particle<T...> &b, const Point &r2) const {
                        return (r2.squaredNorm()<d2(a.id,b.id)) ? pc::infty : 0;
                    }

                void to_json(json &j) const override { j["comment"] = "N/A"; }
                void from_json(const json&) override {}
            }; //!< Hardsphere potential

    }//end of namespace Potential
}//end of namespace Faunus
