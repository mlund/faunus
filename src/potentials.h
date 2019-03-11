#pragma once

#include <array>
#include "geometry.h"
#include "core.h"
#include "auxiliary.h"
#include "tabulate.h"
#include "functionparser.h"

namespace Faunus {
    namespace Potential {

        using namespace std::string_literals;

        struct PairPotentialBase {
            std::string name;
            std::string cite;
            virtual void to_json(json&) const=0;
            virtual void from_json(const json&)=0;
            virtual ~PairPotentialBase();
        }; //!< Base for all pair-potentials

        void to_json(json &j, const PairPotentialBase &base); //!< Serialize any pair potential to json
        void from_json(const json &j, PairPotentialBase &base); //!< Serialize any pair potential from json

        /**
         * @brief Statically combines two pair potentials at compile-time
         *
         * This is the most efficient way to combining pair-potentials due
         * to the possibility for compile-time optimisation.
         */
        template<class T1, class T2>
            struct CombinedPairPotential : public PairPotentialBase {
                T1 first;  //!< First pair potential of type T1
                T2 second; //!< Second pair potential of type T2
                CombinedPairPotential(const std::string &name="") {
                    this->name = name;
                }
                template<typename... T>
                    inline double operator()(const Particle<T...> &a, const Particle<T...> &b, const Point &r) const {
                        return first(a, b, r) + second(a, b, r);
                    } //!< Combine pair energy

                template<typename... T>
                    inline Point force(const Particle<T...> &a, const Particle<T...> &b, double r2, const Point &p) {
                        return first.force(a, b, r2, p) + second.force(a, b, r2, p);
                    } //!< Combine force

                void from_json(const json &j) override {
                    first = j;
                    second = j;
                }

                void to_json(json &j) const override { j = {first,second}; }
            };

        template<class T1, class T2,
            class = typename std::enable_if<std::is_base_of<PairPotentialBase, T1>::value>::type,
            class = typename std::enable_if<std::is_base_of<PairPotentialBase, T2>::value>::type>
                CombinedPairPotential<T1, T2> &operator+(const T1 &pot1, const T2&) {
                    return *(new CombinedPairPotential<T1, T2>(pot1.name));
                } //!< Statically add two pair potentials at compile-time

        struct Dummy : public PairPotentialBase {
            Dummy();
            template<typename... T>
                double operator()(const Particle<T...>&, const Particle<T...>&, const Point&) const {
                    return 0;
                }
            void from_json(const json&) override {}
            void to_json(json&) const override {}
        }; //!< A dummy pair potential that always returns zero

        template<typename Tparticle>
            struct SigmaEpsilonTable {
                enum Mixers {LB, NONE};
                Mixers mixer = NONE;
                PairMatrix<double> s2,eps; // matrix of sigma_ij^2 and 4*eps_ij
            }; //!< Table of sigma and epsilons

        template<typename Tparticle>
            void from_json(const json &j, SigmaEpsilonTable<Tparticle> &m) {
                std::function<std::pair<double,double>(double,double,double,double)> mixerFunc;

                auto mixer = j.at("mixing").get<std::string>();
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
                size_t n=atoms.size(); // number of atom types
                m.s2.resize(n); // not required...
                m.eps.resize(n);// ...but possible reduced mem. fragmentation
                for (auto &i : atoms)
                    for (auto &j : atoms) {
                        double sigma, epsilon; // mixed values
                        std::tie( sigma, epsilon ) = mixerFunc(i.sigma, j.sigma, i.eps, j.eps);
                        m.s2.set(  i.id(), j.id(), sigma*sigma );
                        m.eps.set( i.id(), j.id(), 4*epsilon ); // should already be in kT
                    }

                // custom eps/sigma for specific pairs
                if (j.count("custom")==1) {
                    auto &_j = j.at("custom");
                    if (_j.is_object()) {
                        for (auto it=_j.begin(); it!=_j.end(); ++it) {
                            auto v = words2vec<std::string>( it.key() );
                            if (v.size()==2) {
                                int id1 = (*findName( atoms, v[0])).id();
                                int id2 = (*findName( atoms, v[1])).id();
                                m.s2.set( id1, id2, std::pow( it.value().at("sigma").get<double>(), 2) );
                                m.eps.set(id1, id2, 4*it.value().at("eps").get<double>() * 1.0_kJmol);
                            } else
                                throw std::runtime_error("custom epsilon/sigma parameters require exactly two space-separated atoms");
                        }
                    } else
                        throw std::runtime_error("custom sigma/epsilon syntax error");
                }
            }

        template<typename Tparticle>
            void to_json(json &j, const SigmaEpsilonTable<Tparticle> &m) {
                j["mixing"] = "LB";
                j["epsilon unit"] = "kJ/mol";
                auto& _j = j["custom"];
                for (size_t i=0; i<m.eps.size(); i++)
                    for (size_t j=0; j<m.eps.size(); j++)
                        if (i>=j) {
                            auto str = atoms[i].name+" "+atoms[j].name;
                            _j[str] = { {"eps", m.eps(i,j)/4.0_kJmol}, {"sigma", std::sqrt(m.s2(i,j))}  };
                            _roundjson(_j[str], 5);
                        }
            }

        /**
         * @brief Lennard-Jones with arbitrary mixing rule
         * @note Mixing data is _shared_ upon copying
         */
        template<typename Tparticle>
            class LennardJones : public PairPotentialBase {
                protected:
                    std::shared_ptr<SigmaEpsilonTable<Tparticle>> m; // table w. sigma_ij^2 and 4xepsilon

                public:
                    LennardJones(const std::string &name="lennardjones"s) {
                        PairPotentialBase::name=name;
                        m = std::make_shared<SigmaEpsilonTable<Tparticle>>();
                    }

                    template<typename... T>
                        Point force(const Particle<T...> &a, const Particle<T...> &b, double r2, const Point &p) const {
                            double s6=powi<3>( m->s2(a.id,b.id) );
                            double r6=r2*r2*r2;
                            double r14=r6*r6*r2;
                            return 6.*m->eps(a.id,b.id) * s6 * (2*s6-r6) / r14 * p;
                        }

                    template<typename... T>
                        double operator()(const Particle<T...> &a, const Particle<T...> &b, const Point &r) const {
                            double x=m->s2(a.id,b.id)/r.squaredNorm(); //s2/r2
                            x=x*x*x; // s6/r6
                            return m->eps(a.id,b.id) * (x*x - x);
                        }

                    void to_json(json &j) const override { j = *m; }
                    void from_json(const json &j) override { *m = j; }
            };

        /**
         * @brief Weeks-Chandler-Andersen pair potential
         * @details This is a Lennard-Jones type potential, cut and shifted to zero
         * at @f$r_c=2^{1/6}\sigma@f$. More info can be found in at
         * <http://doi.org/ct4kh9> and the functional form is:
         * @f[
         * \beta u = 4 \beta \epsilon \left ( (b/r)^{12} - (b/r)^6 + \frac{1}{4} \right )
         * @f]
         * where sigma, epsilon per default are set
         * using Lorentz-Berthelot mixing rules.
         *
         * @note Mixing data is _shared_ upon copying
         */
        template<class Tparticle>
            class WeeksChandlerAndersen : public LennardJones<Tparticle> {
                private:
                    typedef LennardJones<Tparticle> base;
                    using base::m;
                    static constexpr double onefourth=0.25, twototwosixth=1.2599210498948732;
                public:
                    WeeksChandlerAndersen(const std::string &name="wca") {
                        base::name=name;
                        base::cite="doi:ct4kh9";
                    }

                    template<typename... T>
                        inline double operator() (const Particle<T...> &a, const Particle<T...> &b, double r2) const {
                            double x=m->s2(a.id,b.id); // s^2
                            if (r2>x*twototwosixth)
                                return 0;
                            x=x/r2;  // (s/r)^2
                            x=x*x*x;// (s/r)^6
                            return m->eps(a.id,b.id)*(x*x - x + onefourth);
                        }

                    template<typename... T>
                        double operator()(const Particle<T...> &a, const Particle<T...> &b, const Point &r) const {
                            return operator()(a,b,r.squaredNorm());
                        }

                    template<typename... T>
                        Point force(const Particle<T...> &a, const Particle<T...> &b, double r2, const Point &p) const {
                            double x=m->s2(a.id,b.id); // s^2
                            if (r2>x*twototwosixth)
                                return Point(0,0,0);
                            x=x/r2;  // (s/r)^2
                            x=x*x*x;// (s/r)^6
                            return m->eps(a.id,b.id)*6*(2*x*x - x)/r2*p;
                        }
            }; // Weeks-Chandler-Andersen potential

        /**
         * @brief Pairwise SASA potential calculating the surface area of inter-secting spheres 
         */
        class SASApotential : public PairPotentialBase {
            private:
                bool shift=true; // shift potential to zero at large separations?
                double proberadius=0, conc=0;

                double area(double R, double r, double d_squared) const; //!< Total surface area of two intersecting spheres or radii R and r as a function of separation

            public:
                template<class Tparticle>
                    double operator()(const Tparticle &a, const Tparticle &b, const Point &r_ab) const {
                        double tfe = 0.5 * ( atoms[a.id].tfe + atoms[b.id].tfe );
                        double tension = 0.5 * ( atoms[a.id].tension + atoms[b.id].tension );
                        if (fabs(tfe)>1e-6 or fabs(tension)>1e-6)
                            return (tension + conc*tfe) *  area(
                                    0.5*atoms[a.id].sigma,
                                    0.5*atoms[b.id].sigma, r_ab.squaredNorm() );
                        return 0;
                    }

                SASApotential(const std::string &name="sasa");
                void to_json(json &j) const override;
                void from_json(const json &j) override;
        };

        struct Coulomb : public PairPotentialBase {
            Coulomb(const std::string &name="coulomb");
            double lB; //!< Bjerrum length
            template<typename... T>
                double operator()(const Particle<T...> &a, const Particle<T...> &b, const Point &r) const {
                    return lB * a.charge * b.charge / r.norm();
                }
            void to_json(json &j) const override;
            void from_json(const json &j) override;
        };

        /**
         * @brief Hardsphere potential
         * @note `PairMatrix` is _shared_ upon copying
         */
        template<class Tparticle>
            class HardSphere : public PairPotentialBase {
                private:
                    std::shared_ptr<PairMatrix<double>> d2; // matrix of (r1+r2)^2
                public:
                    HardSphere(const std::string &name="hardsphere") {
                        PairPotentialBase::name=name;
                        d2 = std::make_shared<PairMatrix<double>>();
                        for (auto &i : atoms)
                            for (auto &j : atoms)
                                d2->set( i.id(), j.id(), std::pow((i.sigma+j.sigma)/2,2));
                    }
                    double operator()(const Tparticle &a, const Tparticle &b, const Point &r) const {
                        return r.squaredNorm() < d2->operator()(a.id,b.id) ? pc::infty : 0;
                    }
                    void to_json(json&) const override {}
                    void from_json(const json&) override {}
            }; //!< Hardsphere potential

        struct RepulsionR3 : public PairPotentialBase {
            double f=0, s=0, e=0;

            RepulsionR3(const std::string &name="repulsionr3");
            void from_json(const json &j) override;
            void to_json(json &j) const override;

            template<class Tparticle>
                double operator() (const Tparticle&, const Tparticle&, const Point &_r) const {
                    double r2 = _r.squaredNorm(), r = sqrt(r2);
                    return f / (r*r2) + e * std::pow( s/r, 12 );
                }
        };

        /**
         * @brief Cosine attraction
         * @details This is an attractive potential used for coarse grained lipids
         * and has the form:
         * @f[
         *     \beta u(r) = -\epsilon \cos^2 [ \pi(r-r_c)/2w_c ]
         * @f]
         * for \f$r_c\leq r \leq r_c+w_c\f$. For \f$r<r_c\f$, \f$\beta u=-\epsilon\f$,
         * while zero for \f$r>r_c+w_c\f$.
         *
         * JSON keywords:
         *
         * Key     | Description
         * :-------| :---------------------------
         * `eps`   | Depth, \f$\epsilon\f$ [kJ/mol]
         * `rc`    | Width, r_c [angstrom]
         * `wc`    | Decay range, w_c [angstrom]
         *
         */
        class CosAttract : public PairPotentialBase {
            double eps, wc, rc, rc2, c, rcwc2;
            public:
            CosAttract(const std::string &name="cos2");

            /**
             * @todo
             * The function `x(c,r2,rc)` defined herein could be approximated
             * by a series expansion for `r2=rcwc2`. In this way one can
             * avoid `cos()` and `sqrt()`. C code for this could be generated
             * in Matlab:
             *
             * ~~~
             * with(CodeGeneration)
             * x := series(cos(c*(sqrt(r2)-rc)), r2 = rcwc2, 2)
             * convert(x, polynom)
             * C(%, resultname = "x")
             * ~~~
             */
            template<typename... T>
                double operator()(const Particle<T...>&, const Particle<T...>&, const Point &r) const {
                    double r2=r.squaredNorm();
                    if (r2<rc2)
                        return -eps;
                    if (r2>rcwc2)
                        return 0;
                    double x=std::cos( c*( sqrt(r2)-rc ) );
                    return -eps*x*x;
                }

            template<class Tparticle>
                Point force(const Tparticle&, const Tparticle&, double r2, const Point &p) {
                    if (r2<rc2 || r2>rcwc2)
                        return Point(0,0,0);
                    double r=sqrt(r2);
                    double x1=std::cos( c*( r-rc ) );
                    double x2=std::sin( c*( r-rc ) );
                    return -2*c*eps*x1*x2/r*p;
                }

            void to_json(json &j) const override;
            void from_json(const json &j) override;
        };

        /**
         * @brief Charge-nonpolar pair interaction
         * @note Pair data is _shared_ upon copying
         */
        template<class Tparticle>
            class Polarizability : public Coulomb {
                private:
                    double epsr;
                    std::shared_ptr<PairMatrix<double>> m_neutral, m_charged;

                public:
                    Polarizability (const std::string &name="polar") {
                        PairPotentialBase::name=name;
                        m_neutral = std::make_shared<PairMatrix<double>>();
                        m_charged = std::make_shared<PairMatrix<double>>();
                    }

                    void from_json(const json &j) override {
                        epsr = j.at("epsr").get<double>();
                        double lB = pc::lB(epsr);
                        for (auto &i : atoms) {
                            for (auto &j : atoms) {
                                m_neutral->set(i.id(), j.id(), -3*i.alphax*pow(0.5*i.sigma,3)*j.alphax*pow(0.5*j.sigma,3) );
                                // titrating particles must be charged in the beginning
                                m_charged->set(i.id(), j.id(), -lB/2 * ( pow(i.charge,2)*j.alphax*pow(0.5*j.sigma,3) +
                                            pow(j.charge,2)*i.alphax*pow(0.5*i.sigma,3) ) );
                            }
                        }
                    }

                    void to_json(json &j) const override {
                        j = { {"epsr",epsr} };
                    }

                    double operator() (const Tparticle &a, const Tparticle &b, const Point &r) const {
                        double r2=r.squaredNorm();
                        double r4inv=1/(r2*r2);
                        if (fabs(a.charge)>1e-9 or fabs(b.charge)>1e-9)
                            return (*m_charged)(a.id,b.id)*r4inv;
                        else
                            return (*m_neutral)(a.id,b.id)/r2*r4inv;
                    }

                    Point force(const Tparticle &a, const Tparticle &b, double r2, const Point &p) {
                        double r6inv=1/(r2*r2*r2);
                        if (fabs(a.charge)>1e-9 or fabs(b.charge)>1e-9)
                            return 4*m_charged->operator()(a.id,b.id)*r6inv*p;
                        else
                            return 6*m_neutral->operator()(a.id,b.id)/r2*r6inv*p;
                    }
            };

        /**
         * @brief Finite Extensible Nonlinear Elastic (FENE) potential
         *
         * This is an anharmonic bonding potential with the form:
         * @f[
         *     \beta u(r) = -\frac{k r_0^2}{2}\ln \left [ 1-(r/r_0)^2 \right ]
         * @f]
         * for \f$r<r_0\f$, otherwise infinity. JSON keywords:
         *
         * - `stiffness` Bond stiffness, `k` [kT]
         * - `maxsep` Maximum separation, `r_0` [angstrom]
         *
         * More info: doi:10.1103/PhysRevE.59.4248
         */
        class FENE : public PairPotentialBase {
            double k, r02, r02inv;

            public:
            FENE(const std::string &name="fene");
            void from_json(const json &j) override;
            void to_json(json &j) const override;

            /** @brief Energy in kT between two particles, r2 = squared distance */
            template<class Tparticle>
                inline double operator() (const Tparticle&, const Tparticle&, const Point &r) {
                    double r2=r.squaredNorm();
                    return (r2>r02) ? pc::infty : -0.5*k*r02*std::log(1-r2*r02inv);
                }

            template<class Tparticle>
                Point force(const Tparticle&, const Tparticle&, double r2, const Point &p) {
                    return (r2>r02) ? -pc::infty*p : -k * r02 / (r02-r2) * p;
                }
        };

        /** @brief Coulomb type potentials with spherical cutoff */
        class CoulombGalore : public PairPotentialBase {
            Tabulate::Andrea<double> sf; // splitting function
            Tabulate::TabulatorBase<double>::data table; // data for splitting function
            std::function<double(double)> calcDielectric; // function for dielectric const. calc.
            std::string type;
            double selfenergy_prefactor;
            double lB, depsdt, rc, rc2, rc1i, epsr, epsrf, alpha, kappa, I;
            int order;
	    bool ellipse_cutoff;
	    Point Rc_ellipse;

            void sfYukawa(const json &j);
            void sfReactionField(const json &j);
            void sfQpotential(const json &j);
            void sfYonezawa(const json &j);
            void sfFanourgakis(const json &j);
            void sfFennel(const json &j);
            void sfEwald(const json &j);
            void sfWolf(const json &j);
            void sfPlain(const json &j, double val=1);

            public:
            CoulombGalore(const std::string &name="coulomb");

            void from_json(const json &j) override;

            template<class Tparticle>
                double operator()(const Tparticle &a, const Tparticle &b, double r2) const {
                    if (r2 < rc2) {
                        double r = std::sqrt(r2);
                        return lB * a.charge * b.charge / r * sf.eval( table, r*rc1i );
                    }
                    return 0;
                }

            template<typename... T>
                double operator()(const Particle<T...> &a, const Particle<T...> &b, const Point &r) const {
                    return operator()(a,b,r.squaredNorm());
                }

            template<typename... T>
                Point force(const Particle<T...> &a, const Particle<T...> &b, double r2, const Point &p) const {
                    if (r2 < rc2) {
                        double r = sqrt(r2);
                        return lB * a.charge * b.charge * ( -sf.eval( table, r*rc1i )/r2 + sf.evalDer( table, r*rc1i )/r )*p;
                    }
                    return Point(0,0,0);
                }

            /**
             * @brief Self-energy of the potential
             */
            template<class Tpvec, class Tgroup>
                double internal(const Tgroup &g) const {
                    double Eq = 0;
                    for (auto i : g)
                        Eq += i.charge * i.charge;
                    return -selfenergy_prefactor*Eq*lB/rc;
                }

            double dielectric_constant(double M2V);

            void to_json(json &j) const override;
        };

        /**
         * @brief Custom pair-potential taking math. expressions at runtime
         */
        class CustomPairPotential : public PairPotentialBase {
            private:
                ExprFunction<double> expr;
                struct Data {
                    double r=0, q1=0, q2=0, s1=0, s2=0;
                };
                double Rc2;
                std::shared_ptr<Data> d;
                json jin; // initial json input
            public:
                template<typename... T>
                    inline double operator()(const Particle<T...> &a, const Particle<T...> &b, const Point &r) const {
                        double r2 = r.squaredNorm();
                        if (r2>Rc2)
                            return 0;
                        d->r  = sqrt(r2);
                        d->q1 = a.charge;
                        d->q2 = b.charge;
                        d->s1 = atoms[a.id].sigma;
                        d->s2 = atoms[b.id].sigma;
                        return expr();
                    }
                CustomPairPotential(const std::string &name="custom");
                void from_json(const json&) override;
                void to_json(json&) const override;
        };

#ifdef DOCTEST_LIBRARY_INCLUDED
        TEST_CASE("[Faunus] CustomPairPotential")
        {
            using doctest::Approx;
            typedef Particle<Radius, Charge, Dipole, Cigar> Tparticle;
            json j = R"({ "atomlist" : [
                 {"A": { "q":1.0,  "r":3, "eps":0.1 }},
                 {"B": { "q":-1.0, "r":4, "eps":0.05 }} ]})"_json;

            atoms = j["atomlist"].get<decltype(atoms)>();

            Tparticle a,b;
            a = atoms[0];
            b = atoms[1];

            CustomPairPotential pot = R"({
                "constants": { "kappa": 30, "lB": 7},
                "function": "lB * q1 * q2 / (s1+s2) * exp(-kappa/r) * kT + pi"})"_json;

            CHECK( pot(a,b,{0,0,2}) == Approx( -7/(3.0+4.0) * std::exp(-30/2) * pc::kT() + pc::pi  ));
        }
#endif

        /**
         * @brief Arbitrary potentials for specific atom types
         *
         * This maintains a species x species matrix with function pointers (`std::function`)
         * that wraps pair potentials. Flexibility over performance.
         *
         * @todo `to_json` should retrieve info from potentials instead of merely passing input
         * @warning Each atom pair will be assigned an instance of a pair-potential. This *could* be
         *          problematic if these have large memory requirements.
         */
        template<class T /** particle type */>
            class FunctorPotential : public PairPotentialBase {
                typedef std::function<double(const T&, const T&, const Point&)> uFunc;
                PairMatrix<uFunc,true> umatrix; // matrix with potential for each atom pair
                typedef Tabulate::TabulatorBase<double>::data Ttable; // data for tabulated potential
                Tabulate::Andrea<double> tblt; // tabulated potential
                PairMatrix<Ttable,true> tmatrix; // matrix with tabulated potential for each atom pair
                json _j; // storage for input json
                double rc2;
                typedef CombinedPairPotential<Coulomb,HardSphere<T>> PrimitiveModel;
                typedef CombinedPairPotential<Coulomb,WeeksChandlerAndersen<T>> PrimitiveModelWCA;

                // List of pair-potential instances used when constructing functors.
                // Note that potentials w. large memory requirements (LJ, WCA etc.)
                // typically use `shared_ptr` so that the created functors _share_
                // the data. That is *only* put the pair-potential here if you can
                // share internal (shared) pointers.
                std::tuple<
                    CoulombGalore,    // 0
                    CosAttract,       // 1
                    Polarizability<T>,// 2
                    HardSphere<T>,    // 3
                    LennardJones<T>,  // 4
                    RepulsionR3,      // 5
                    SASApotential,    // 6
                    WeeksChandlerAndersen<T>,// 7
                    PrimitiveModel,   // 8
                    PrimitiveModelWCA // 9
                        > potlist;

                 uFunc combineFunc(const json &j) { 
                    uFunc u = [](const T&, const T&, const Point&){return 0.0;};
                    if (j.is_array()) {
                        for (auto &i : j) // loop over all defined potentials in array
                            if (i.is_object() and (i.size()==1))
                                for (auto it : i.items()) {
                                    uFunc _u = nullptr;
                                    try {
                                        if (it.key()=="custom") _u = CustomPairPotential() = it.value();
                                        else if (it.key()=="coulomb") _u = std::get<0>(potlist) = i;
                                        else if (it.key()=="cos2") _u = std::get<1>(potlist) = i;
                                        else if (it.key()=="polar") _u = std::get<2>(potlist) = i;
                                        else if (it.key()=="hardsphere") _u = std::get<3>(potlist) = i;
                                        else if (it.key()=="lennardjones") _u = std::get<4>(potlist) = i;
                                        else if (it.key()=="repulsionr3") _u = std::get<5>(potlist) = i;
                                        else if (it.key()=="sasa") _u = std::get<6>(potlist) = i;
                                        else if (it.key()=="wca") _u = std::get<7>(potlist) = i;
                                        else if (it.key()=="pm") _u = std::get<8>(potlist) = it.value();
                                        else if (it.key()=="pmwca") _u = std::get<9>(potlist) = it.value();
                                        // place additional potentials here...
                                    } catch (std::exception &e) {
                                        throw std::runtime_error("Error adding energy '" + it.key() + "': " + e.what() + usageTip[it.key()]);
                                    }

                                    if (_u!=nullptr) // if found, sum them into new function object
                                        u = [u,_u](const T&a, const T&b, const Point &r){return u(a,b,r)+_u(a,b,r);};
                                    else
                                        throw std::runtime_error("unknown pair-potential: " + it.key());
                                }
                    } else
                        throw std::runtime_error("dictionary of potentials required");
                    return u;
                } // parse json array of potentials to a single potential function object

                public:

                FunctorPotential(const std::string &name="") {
                    PairPotentialBase::name = name;
                }

                double operator()(const T &a, const T &b, const Point &r) const {
                    double r2 = r.squaredNorm();
                    if (r2 > tmatrix(a.id, b.id).rmax2)
                        return 0.0;
                    else if (r2 <= tmatrix(a.id, b.id).rmin2)
                        return umatrix(a.id, b.id)(a, b, Point(0,0,sqrt(r2))); // pc::infty;
                    else 
                        return tblt.eval(tmatrix(a.id, b.id), r2);
                }

                void to_json(json &j) const override { j = _j; }

                void from_json(const json &j) override {
                    tblt.setTolerance(j.value("utol",1e-5),j.value("ftol",1e-2) );
                    _j = j;
                    double rmax2 = pc::Nav;
                    umatrix = decltype(umatrix)( atoms.size(), combineFunc(j.at("default")) );
                    for (auto it=j.begin(); it!=j.end(); ++it) {
                        auto atompair = words2vec<std::string>(it.key()); // is this for a pair of atoms?
                        if (atompair.size()==2) {
                            auto ids = names2ids(atoms, atompair);
                            umatrix.set(ids[0], ids[1], combineFunc(it.value()));
                        }
                    }
                    for (size_t i=0; i<atoms.size(); ++i) {
                        for (size_t k=0; k<=i; ++k) {
                           if (atoms[i].implicit==false and atoms[k].implicit==false) {
                                T a = atoms.at(i);
                                T b = atoms.at(k);
                                double rmin2 = .5*(atoms[i].sigma + atoms[k].sigma);
                                rmin2 = rmin2*rmin2;
                                auto it = j.find("cutoff_g2g");
                                if (it->is_number())
                                    rmax2 = std::pow( it->get<double>(), 2 );
                                else if (it->is_object())
                                    rmax2 = std::pow( it->at("default").get<double>(), 2);
                                while (rmin2 >= 1e-2) {
                                    if (std::fabs(umatrix(i,k)(a, b, Point(0,0,sqrt(rmin2)))) > 1e6)
                                        rmin2 = rmin2 + 1e-2;
                                    else if (std::fabs(umatrix(i,k)(a, b, Point(0,0,sqrt(rmin2)))) > 1e5)
                                        break;
                                    else
                                        rmin2 = rmin2 - 1e-2;
                                }
                                while (rmax2 >= 1e-2) {
                                    if (std::fabs(umatrix(i,k)(a, b, Point(0,0,sqrt(rmax2)))) > pc::epsilon_dbl)
                                        break;
                                    rmax2 = rmax2 - 1e-2;
                                }
                                Ttable knotdata = tblt.generate( [&](double r2) { return umatrix(i,k)(a, b, Point(0,0,sqrt(r2))); }, rmin2, rmax2);
                                tmatrix.set(i, k, knotdata);
                                std::ofstream file(atoms[i].name+"-"+atoms[k].name+"_tabulated.dat"); // output file
                                file << "# Separation\tTabulated\tOriginal\n";
                                double r2 = rmin2;
                                while (r2 < rmax2) {
                                    r2 = r2 + 1e-2;
                                    file << sqrt(r2) << "\t" << tblt.eval(tmatrix(i, k), r2) << "\t" << umatrix(i,k)(a, b, Point(0,0,sqrt(r2))) << "\n";
                                }
                            }
                        }
                    }
                }
            };

#ifdef DOCTEST_LIBRARY_INCLUDED
        TEST_CASE("[Faunus] FunctorPotential")
        {
            using doctest::Approx;
            typedef Particle<Radius, Charge, Dipole, Cigar> T;

            json j = R"({ "atomlist" : [
                 {"A": { "q":1.0,  "r":1.1, "eps":0.1 }},
                 {"B": { "q":-1.0, "r":2.0, "eps":0.05 }},
                 {"C": { "r":1.0 }} ]})"_json;

            atoms = j["atomlist"].get<decltype(atoms)>();

            FunctorPotential<T> u = R"(
                {
                  "default": [
                    { "coulomb" : {"epsr": 80.0, "type": "plain", "cutoff":20} }
                  ],
                  "A B" : [
                    { "coulomb" : {"epsr": 80.0, "type": "plain", "cutoff":20} },
                    { "wca" : {"mixing": "LB"} }
                  ],
                  "C C" : [
                    { "hardsphere" : {} }
                  ]
                 }
                )"_json;

            Coulomb coulomb = R"({ "coulomb": {"epsr": 80.0, "type": "plain", "cutoff":20} } )"_json;
            WeeksChandlerAndersen<T> wca = R"({ "wca" : {"mixing": "LB"} })"_json;

            T a = atoms[0];
            T b = atoms[1];
            T c = atoms[2];
            Point r={2,0,0};
            CHECK( u(a,a,r) == Approx( coulomb(a,a,r) ) );
            CHECK( u(b,b,r) == Approx( coulomb(b,b,r) ) );
            CHECK( u(a,b,r) == Approx( coulomb(a,b,r) + wca(a,b,r) ) );
            CHECK( u(c,c,r*1.01) == 0 );
            CHECK( u(c,c,r*0.99) == pc::infty );
        }
#endif

        /**
         * @brief Base class for bonded potentials
         *
         * This stores data on the bond type; atom indices; json keywords;
         * and potentially also the energy function (nullptr per default).
         */
        struct BondData {
            enum Variant {HARMONIC=0, FENE, HARMONIC_TORSION, G96_TORSION, PERIODIC_DIHEDRAL, NONE};
            std::vector<int> index;
            bool exclude=false;           //!< True if exclusion of non-bonded interaction should be attempted 
            bool keepelectrostatics=true; //!< If `exclude==true`, try to keep electrostatic interactions
            std::function<double(Geometry::DistanceFunction)> energy=nullptr; //!< potential energy (kT)

            virtual void from_json(const json&)=0;
            virtual void to_json(json&) const=0;
            virtual int numindex() const=0; //!< Required number of atom indices for bond
            virtual Variant type() const=0; //!< Returns bond type (sett `Variant` enum)
            virtual std::string name() const=0; //!< Name/key of bond type used in for json I/O
            virtual std::shared_ptr<BondData> clone() const=0; //!< Make shared pointer *copy* of data
            bool hasEnergyFunction() const; //!< test if energy function has been set
            void shift( int offset ); //!< Shift indices
            virtual ~BondData();
        };

        /**
         * @brief Harmonic Bond
         */
        struct HarmonicBond : public BondData {
            double k=0, req=0;
            int numindex() const override;
            Variant type() const override;
            std::shared_ptr<BondData> clone() const override;
            void from_json(const json &j) override;
            void to_json(json &j) const override;
            std::string name() const override;

            template<typename Tpvec>
                void setEnergyFunction(const Tpvec &p) {
                    energy = [&](Geometry::DistanceFunction dist) {
                        double d = req - dist(p[index[0]].pos, p[index[1]].pos).norm();
                        return k*d*d;
                    };
                }
        };

        /**
         * @brief FENE bond
         */
        struct FENEBond : public BondData {
            std::array<double,4> k = {{0,0,0,0}};
            int numindex() const override;
            Variant type() const override;
            std::shared_ptr<BondData> clone() const override;
            void from_json(const json &j) override;
            void to_json(json &j) const override;
            std::string name() const override;

            template<typename Tpvec>
                void setEnergyFunction(const Tpvec &p) {
                    energy = [&](Geometry::DistanceFunction dist) {
                        double wca=0, d=dist( p[index[0]].pos, p[index[1]].pos ).squaredNorm();
                        double x = k[3];
                        if (d<=x*1.2599210498948732) {
                            x = x/d;
                            x = x*x*x;
                            wca = k[2]*(x*x - x + 0.25);
                        }
                        return (d>k[1]) ? pc::infty : -0.5*k[0]*k[1]*std::log(1-d/k[1]) + wca;
                    };
                }
        }; // end of FENE

        struct HarmonicTorsion : public BondData {
            double k=0, aeq=0;
            int numindex() const override;
            std::shared_ptr<BondData> clone() const override;
            void from_json(const json &j) override;
            void to_json(json &j) const override;
            Variant type() const override;
            std::string name() const override;

            template<typename Tpvec>
                void setEnergyFunction(const Tpvec &p) {
                    energy = [&](Geometry::DistanceFunction dist) {
                        Point ray1 = dist( p[index[0]].pos, p[index[1]].pos );
                        Point ray2 = dist( p[index[2]].pos, p[index[1]].pos );
                        double angle = std::acos(ray1.dot(ray2)/ray1.norm()/ray2.norm());
                        return 0.5 * k * (angle - aeq) * (angle - aeq);
                    };
                }
        }; // end of HarmonicTorsion

        struct GromosTorsion : public BondData {
            double k=0, aeq=0;
            int numindex() const override;
            std::shared_ptr<BondData> clone() const override;
            void from_json(const json &j) override;
            void to_json(json &j) const override;
            Variant type() const override;
            std::string name() const override;

            template<typename Tpvec>
                void setEnergyFunction(const Tpvec &p) {
                    energy = [&](Geometry::DistanceFunction dist) {
                        Point ray1 = dist( p[index[0]].pos, p[index[1]].pos );
                        Point ray2 = dist( p[index[2]].pos, p[index[1]].pos );
                        double dangle = aeq-std::acos(ray1.dot(ray2)/ray1.norm()/ray2.norm());
                        return k * dangle * dangle;
                    };
                }
        }; // end of GromosTorsion

        struct PeriodicDihedral : public BondData {
            std::array<double,3> k;
            int numindex() const override;
            std::shared_ptr<BondData> clone() const override;
            void from_json(const json &j) override;
            void to_json(json &j) const override;
            Variant type() const override;
            std::string name() const override;

            template<typename Tpvec>
                void setEnergyFunction(const Tpvec &p) {
                    energy = [&](Geometry::DistanceFunction dist) {
                        Point vec1 = dist( p[index[1]].pos, p[index[0]].pos );
                        Point vec2 = dist( p[index[2]].pos, p[index[1]].pos );
                        Point vec3 = dist( p[index[3]].pos, p[index[2]].pos );
                        Point norm1 = vec1.cross(vec2);
                        Point norm2 = vec2.cross(vec3);
                        // atan2( [v1×v2]×[v2×v3]⋅[v2/|v2|], [v1×v2]⋅[v2×v3] )
                        double angle = atan2((norm1.cross(norm2)).dot(vec2)/vec2.norm(), norm1.dot(norm2));
                        return k[0] * (1 + cos(k[1]*angle - k[2]));
                    };
                }
        }; // end of PeriodicDihedral

        /*
         * Serialize to/from json
         */

        void to_json(json &j, const std::shared_ptr<BondData> &b);
        void from_json(const json &j, std::shared_ptr<BondData> &b);

        template<typename Tpvec>
            void setBondEnergyFunction(std::shared_ptr<BondData> &b, const Tpvec &p) {
                if (b->type()==BondData::HARMONIC)
                    std::dynamic_pointer_cast<HarmonicBond>(b)->setEnergyFunction(p);
                else if (b->type()==BondData::FENE)
                    std::dynamic_pointer_cast<FENEBond>(b)->setEnergyFunction(p);
                else {
                    assert(false); // we should never reach here
                }
            } //!< Set the bond energy function of `BondData` which require a reference to the particle vector

        inline auto filterBonds(const std::vector<std::shared_ptr<BondData>> &bonds, BondData::Variant bondtype) {
            std::vector<std::shared_ptr<BondData>> filt;
            filt.reserve(bonds.size());
            std::copy_if(bonds.begin(), bonds.end(), std::back_inserter(filt),
                    [bondtype](const auto &d){return d->type()==bondtype;} );
            return filt;
        } //!< Filter bond container for matching bond type and return _reference_ to original

#ifdef DOCTEST_LIBRARY_INCLUDED
        TEST_CASE("[Faunus] BondData")
        {
            std::shared_ptr<BondData> b;

            // exact match required
            CHECK_THROWS( b = R"({ "harmoNIC": {"index":[2,3], "k":0.5, "req":2.1}} )"_json; );

            // test harmonic
            SUBCASE("HarmonicBond") {
                json j = R"({ "harmonic": {"index":[2,3], "k":0.5, "req":2.1}} )"_json;
                b = j;
                CHECK( j == json(b) );
                CHECK_THROWS( b = R"({"harmonic": { "index":[2], "k":0.5, "req":2.1}} )"_json );
                CHECK_THROWS( b = R"({"harmonic": { "index":[2,3], "req":2.1}} )"_json );
                CHECK_THROWS( b = R"({"harmonic": { "index":[2,3], "k":2.1}} )"_json );
            }

            // test fene
            SUBCASE("FENEBond") {
                json j = R"({"fene": { "index":[2,3], "k":1, "rmax":2.1, "eps":2.48, "sigma":2 }} )"_json;
                b = j;
                CHECK( j == json(b) );
                CHECK_THROWS( b = R"({"fene": { "index":[2,3,4], "k":1, "rmax":2.1}} )"_json );
                CHECK_THROWS( b = R"({"fene": { "index":[2,3], "rmax":2.1}} )"_json );
                CHECK_THROWS( b = R"({"fene": { "index":[2,3], "k":1}} )"_json );
                j = json::object();
                CHECK_THROWS( b = j );
            }

            // test harmonic
            SUBCASE("HarmonicTorsion") {
                json j = R"({ "harmonic_torsion": {"index":[0,1,2], "k":0.5, "aeq":60}} )"_json;
                b = j;
                CHECK( j == json(b) );
                CHECK_THROWS( b = R"({"harmonic_torsion": { "index":[2], "k":0.5, "aeq":2.1}} )"_json );
                CHECK_THROWS( b = R"({"harmonic_torsion": { "index":[0,1,2], "aeq":2.1}} )"_json );
                CHECK_THROWS( b = R"({"harmonic_torsion": { "index":[0,1,3], "k":2.1}} )"_json );
            }

            // test bond filter
            SUBCASE("filterBonds()") {
                std::vector<std::shared_ptr<BondData>> bonds = {
                    R"({"fene":      {"index":[2,3], "k":1, "rmax":2.1, "eps":2.48, "sigma":2 }} )"_json,
                    R"({"harmonic" : {"index":[2,3], "k":0.5, "req":2.1} } )"_json
                };
                auto filt = filterBonds(bonds, BondData::HARMONIC);
                CHECK( filt.size() == 1 );
                CHECK( filt[0]->type() == BondData::HARMONIC);
                CHECK( filt[0] == bonds[1] ); // filt should contain references to bonds
            }
        }

        TEST_CASE("[Faunus] Pair Potentials")
        {
            using doctest::Approx;
            json j = R"({ "atomlist" : [
                 { "A": { "r": 1.5, "tension": 0.023} },
                 { "B": { "r": 2.1, "tfe": 0.98 } }]})"_json;

            typedef Particle<Radius> T;
            atoms = j["atomlist"].get<decltype(atoms)>();
            T a,b;
            a.id = 0;
            b.id = 1;

            SUBCASE("SASApotential") {
                SASApotential pot;
                json in = R"({ "sasa": {"molarity": 1.0, "radius": 0.0, "shift":false}})"_json;
                pot = in["sasa"];
                double conc = 1.0 * 1.0_molar;
                double tension = atoms[a.id].tension / 2;
                double tfe = atoms[b.id].tfe / 2;
                double f = tension + conc*tfe;
                CHECK( tension>0.0 );
                CHECK( conc>0.0 );
                CHECK( tfe>0.0 );
                CHECK( f>0.0 );
                CHECK( in == json(pot) );
                CHECK( pot(a,b,{0,0,0})  == Approx( f * 4*pc::pi*2.1*2.1) ); // complete overlap
                CHECK( pot(a,b,{10,0,0}) == Approx( f * 4*pc::pi*(2.1*2.1+1.5*1.5) ) ); // far apart
                CHECK( pot(a,b,{2.5,0,0})== Approx( f * 71.74894965974514 ) ); // partial overlap
            }
        }
#endif

    }//end of namespace Potential
}//end of namespace Faunus
