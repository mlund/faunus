#pragma once

#include "core.h"
#include "auxiliary.h"
#include "tabulate.h"
#include "multipole.h"

namespace Faunus {
    namespace Potential {

        using namespace std::string_literals;

        struct PairPotentialBase {
            std::string name;
            std::string cite;
            virtual void to_json(json&) const=0;
            virtual void from_json(const json&)=0;
        }; //!< Base for all pair-potentials

        void to_json(json &j, const PairPotentialBase &base) {
            base.name.empty() ? base.to_json(j) : base.to_json(j[base.name]);
        } //!< Serialize any pair potential to json

        void from_json(const json &j, PairPotentialBase &base) {
            if (!base.name.empty())
                if (j.count(base.name)==1) {
                    base.from_json(j.at(base.name));
                    return;
                }
            base.from_json(j);
        } //!< Serialize any pair potential from json

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
                    return *(new CombinedPairPotential<T1, T2>(pot1.name));
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
                if (j.count("custom")==1) {
                    auto &_j = j.at("custom");
                    if (_j.is_object()) {
                        for (auto it=_j.begin(); it!=_j.end(); ++it) {
                            auto v = words2vec<std::string>( it.key() );
                            if (v.size()==2) {
                                int id1 = (*findName( atoms<Tparticle>, v[0])).id();
                                int id2 = (*findName( atoms<Tparticle>, v[1])).id();
                                m.s2.set( id1, id2, std::pow( it.value().at("sigma").get<double>(), 2) );
                                m.eps.set(id1, id2, 4*it.value().at("eps").get<double>() * 1.0_kJmol);
                            } else
                                std::runtime_error("custom epsilon/sigma parameters require exactly two space-separated atoms");
                        }
                    }
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
                            auto str = atoms<Tparticle>[i].name+" "+atoms<Tparticle>[j].name;
                            _j[str] = { {"eps", m.eps(i,j)/4.0_kJmol}, {"sigma", std::sqrt(m.s2(i,j))}  };
                            _roundjson(_j[str], 5);
                        }
            }

        /**
         * @brief Lennard-Jones with arbitrary mixing rule
         */
        template<typename Tparticle>
            struct LennardJones : public PairPotentialBase {
                LennardJones(const std::string &name="lennardjones"s) { PairPotentialBase::name=name; }
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
         */
        template<class Tparticle>
            class WeeksChandlerAndersen : public LennardJones<Tparticle> {
                private:
                    typedef LennardJones<Tparticle> base;
                    using base::m;
                    static constexpr double onefourth=0.25, twototwosixth=1.2599210498948732;
                public:
                    inline WeeksChandlerAndersen(const std::string &name="wca") {
                        base::name=name;
                        base::cite="doi:ct4kh9";
                    }

                    template<typename... T>
                        inline double operator() (const Particle<T...> &a, const Particle<T...> &b, double r2) const {
                            double x=m.s2(a.id,b.id); // s^2
                            if (r2>x*twototwosixth)
                                return 0;
                            x=x/r2;  // (s/r)^2
                            x=x*x*x;// (s/r)^6
                            return m.eps(a.id,b.id)*(x*x - x + onefourth);
                        }

                    template<typename... T>
                        double operator()(const Particle<T...> &a, const Particle<T...> &b, const Point &r) const {
                            return operator()(a,b,r.squaredNorm());
                        }

                    template<typename... T>
                        Point force(const Particle<T...> &a, const Particle<T...> &b, double r2, const Point &p) const {
                            double x=m.s2(a.id,b.id); // s^2
                            if (r2>x*twototwosixth)
                                return Point(0,0,0);
                            x=x/r2;  // (s/r)^2
                            x=x*x*x;// (s/r)^6
                            return m.eps(a.id,b.id)*6*(2*x*x - x)/r2*p;
                        }
            }; // Weeks-Chandler-Andersen potential


        struct Coulomb : public PairPotentialBase {
            Coulomb(const std::string &name="coulomb") { PairPotentialBase::name=name; }
            double lB; //!< Bjerrum length
            template<typename... T>
                double operator()(const Particle<T...> &a, const Particle<T...> &b, const Point &r) const {
                    return lB * a.charge * b.charge / r.norm();
                }
            void to_json(json &j) const override { j["epsr"] = pc::lB(lB); }
            void from_json(const json &j) override { lB = pc::lB( j.at("epsr") ); }
        };

        template<class Tparticle>
            struct HardSphere : public PairPotentialBase {
                PairMatrix<double> d2; // matrix of (r1+r2)^2
                HardSphere(const std::string &name="hardsphere") {
                    PairPotentialBase::name=name;
                    for (auto &i : atoms<Tparticle>)
                        for (auto &j : atoms<Tparticle>)
                            d2.set( i.id(), j.id(), std::pow((i.sigma+j.sigma)/2,2));
                }
                double operator()(const Tparticle &a, const Tparticle &b, const Point &r) const {
                    return r.squaredNorm() < d2(a.id,b.id) ? pc::infty : 0;
                }
                void to_json(json &j) const override {}
                void from_json(const json&) override {}
            }; //!< Hardsphere potential

        struct RepulsionR3 : public PairPotentialBase {
            double f=0, s=0, e=0;
            inline RepulsionR3(const std::string &name="repulsionr3") {
                PairPotentialBase::name = name;
            }
            void from_json(const json &j) override {
                f = j.value("prefactor", 1.0);
                e = j.value("lj-prefactor", 1.0);
                s = j.value("sigma", 1.0);
            }
            void to_json(json &j) const override {
                j = {{"prefactor",f}, {"lj-prefactor", e},{"sigma",s}};
            }
            template<class Tparticle>
                double operator() (const Tparticle &a, const Tparticle &b, const Point &_r) const {
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
            CosAttract(const std::string &name="cos2") { PairPotentialBase::name=name; }

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
                double operator()(const Particle<T...> &a, const Particle<T...> &b, const Point &r) const {
                    double r2=r.squaredNorm();
                    if (r2<rc2)
                        return -eps;
                    if (r2>rcwc2)
                        return 0;
                    double x=std::cos( c*( sqrt(r2)-rc ) );
                    return -eps*x*x;
                }

            template<class Tparticle>
                Point force(const Tparticle &a, const Tparticle &b, double r2, const Point &p) {
                    if (r2<rc2 || r2>rcwc2)
                        return Point(0,0,0);
                    double r=sqrt(r2);
                    double x1=std::cos( c*( r-rc ) );
                    double x2=std::sin( c*( r-rc ) );
                    return -2*c*eps*x1*x2/r*p;
                }

            void to_json(json &j) const override {
                j = {{"eps",eps / 1.0_kJmol}, {"rc",rc / 1.0_angstrom }, {"wc", wc / 1.0_angstrom }};
            }

            void from_json(const json &j) override {
                eps = j.at("eps").get<double>() * 1.0_kJmol;
                rc = j.at("rc").get<double>() * 1.0_angstrom ;
                wc = j.at("wc").get<double>() * 1.0_angstrom ;
                rc2 = rc * rc;
                c = pc::pi / 2 / wc;
                rcwc2 = pow((rc + wc), 2);
            }
        };

        /**
         * @brief Charge-nonpolar pair interaction
         */
        template<class Tparticle>
            class Polarizability : public Coulomb {
                private:
                    double epsr;
                    PairMatrix<double> m_neutral, m_charged; 

                public:
                    Polarizability (const std::string &name="polar") { PairPotentialBase::name=name; }

                    inline void from_json(const json &j) override {
                        epsr = j.at("epsr").get<double>();
                        double lB = pc::lB(epsr);
                        for (auto &i : atoms<Tparticle>) {
                            for (auto &j : atoms<Tparticle>) { 
                                m_neutral.set(i.id(), j.id(), -3*i.alphax*pow(0.5*i.sigma,3)*j.alphax*pow(0.5*j.sigma,3) );
                                // titrating particles must be charged in the beginning
                                m_charged.set(i.id(), j.id(), -lB/2 * ( pow(i.p.charge,2)*j.alphax*pow(0.5*j.sigma,3) +
                                        pow(j.p.charge,2)*i.alphax*pow(0.5*i.sigma,3) ) );
                            }
                        }
                    }

                    inline void to_json(json &j) const override {
                        j = { {"epsr",epsr} };
                    }

                    double operator() (const Tparticle &a, const Tparticle &b, const Point &r) const {
                        double r2=r.squaredNorm();
                        double r4inv=1/(r2*r2);
                        if (fabs(a.charge)>1e-9 or fabs(b.charge)>1e-9) 
                            return m_charged(a.id,b.id)*r4inv;
                        else
                            return m_neutral(a.id,b.id)/r2*r4inv;
                    }

                    Point force(const Tparticle &a, const Tparticle &b, double r2, const Point &p) {
                        double r6inv=1/(r2*r2*r2);
                        if (fabs(a.charge)>1e-9 or fabs(b.charge)>1e-9) 
                            return 4*m_charged(a.id,b.id)*r6inv*p;
                        else
                            return 6*m_neutral(a.id,b.id)/r2*r6inv*p;
                    }
            };

        template<class Tparticle>
            class DesernoMembrane : public PairPotentialBase {

                WeeksChandlerAndersen<Tparticle> wca;
                CosAttract cos2;
                int tail;

                public:
                DesernoMembrane(const std::string &name="dmembrane") {
                    PairPotentialBase::name=name;
                    PairPotentialBase::name.clear();
                }

                inline void from_json(const json &j) override {
                    wca = j;
                    cos2 = j;
                    auto it = findName(atoms<Tparticle>, "TL");
                    if ( it!=atoms<Tparticle>.end() )
                        tail = it->id();
                    else
                        throw std::runtime_error("Atom type 'TL' is not defined.");
                }
                void to_json(json &j) const override {
                    json _j;
                    wca.to_json(j);
                    cos2.to_json(_j);
                    j = merge(j,_j);
                }

                double operator() (const Tparticle &a, const Tparticle &b, const Point &r) const {
                    double u=wca(a,b,r);
                    if (a.id==tail and b.id==tail) 
                        u+=cos2(a,b,r);
                    return u;
                }
            };

        template<class Tparticle>
            class DesernoMembraneAA : public PairPotentialBase {

                WeeksChandlerAndersen<Tparticle> wca;
                CosAttract cos2;
                Polarizability<Tparticle> polar;
                int tail;
                int aa;

                public:
                DesernoMembraneAA(const std::string &name="dmembraneAA") {
                    PairPotentialBase::name=name;
                    PairPotentialBase::name.clear();
                }

                inline void from_json(const json &j) override {
                    wca = j;
                    cos2 = j;
                    polar = j;
                    auto it = findName(atoms<Tparticle>, "TL");
                    if ( it!=atoms<Tparticle>.end() )
                        tail = it->id();
                    else
                        throw std::runtime_error("Atom type 'TL' is not defined.");
                    it = findName(atoms<Tparticle>, "AA");
                    if ( it!=atoms<Tparticle>.end() )
                        aa = it->id();
                    else
                        throw std::runtime_error("Atom type 'AA' is not defined.");

                }
                void to_json(json &j) const override {
                    json _j;
                    wca.to_json(j);
                    cos2.to_json(_j);
                    j = merge(j,_j);
                    polar.to_json(_j);
                    j = merge(j,_j);
                }

                double operator() (const Tparticle &a, const Tparticle &b, const Point &r) const {
                    double u=wca(a,b,r);
                    if (a.id==tail and b.id==tail)
                        u+=cos2(a,b,r);
                    if (a.id==aa or b.id==aa) {
                        u+=polar(a,b,r);
                    }
                    return u;
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
            FENE(const std::string &name="fene") { PairPotentialBase::name=name; }

            inline void from_json(const json &j) override {
                k  = j.at("stiffness");
                r02 = std::pow( double(j.at("maxsep")), 2);
                r02inv = 1/r02;
            }

            inline void to_json(json &j) const override {
                j = {{"stiffness",k}, {"maxsep",std::sqrt(r02)}};
            }

            /** @brief Energy in kT between two particles, r2 = squared distance */
            template<class Tparticle>
                inline double operator() (const Tparticle &a, const Tparticle &b, const Point &r) {
                    double r2=r.squaredNorm();
                    return (r2>r02) ? pc::infty : -0.5*k*r02*std::log(1-r2*r02inv);
                }

            template<class Tparticle>
                Point force(const Tparticle &a, const Tparticle &b, double r2, const Point &p) {
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

            void sfYukawa(const json &j) {
                kappa = 1.0 / j.at("debyelength").get<double>();
                I = kappa*kappa / ( 8.0*lB*pc::pi*pc::Nav/1e27 );
                table = sf.generate( [&](double q) { return std::exp(-q*rc*kappa) - std::exp(-kappa*rc); }, 0, 1 ); // q=r/Rc 
                // we could also fill in some info std::string or JSON output...
            }

            void sfReactionField(const json &j) {
                epsrf = j.at("eps_rf");
                table = sf.generate( [&](double q) { return 1 + (( epsrf - epsr ) / ( 2 * epsrf + epsr ))*q*q*q
                        - 3 * ( epsrf / ( 2 * epsrf + epsr ))*q ; }, 0, 1); 
                calcDielectric = [&](double M2V) {
                    if(epsrf > 1e10)
                        return 1 + 3*M2V;
                    if(fabs(epsrf-epsr) < 1e-6)
                        return 2.25*M2V + 0.25 + 0.75*sqrt(9*M2V*M2V + 2*M2V + 1);
                    if(fabs(epsrf-1.0) < 1e-6)
                        return ( 2*M2V + 1 ) / ( 1 - M2V );
                    return 0.5 * ( 2*epsrf - 1 + sqrt( -72*M2V*M2V*epsrf
                                + 4*epsrf*epsrf + 4*epsrf + 1) ) / ( 3*M2V-1 ); // Needs to be checked!
                    //return (6*M2V*epsrf + 2*epsrf + 1.0)/(1.0 + 2*epsrf - 3*M2V); // Is OK when epsr=1.0
                };
                selfenergy_prefactor = 1.5*epsrf/(2.0*epsrf + epsr); // Correct?!, see Eq.14 in DOI: 10.1021/jp510612w
                // we could also fill in some info std::string or JSON output...
            }

            void sfQpotential(const json &j)
            {
                order = j.value("order",300);
                table = sf.generate( [&](double q) { return qPochhammerSymbol( q, 1, order ); }, 0, 1 );
                calcDielectric = [&](double M2V) { return 1 + 3*M2V; };
                selfenergy_prefactor = 0.5;
            }

            void sfYonezawa(const json &j)
            {
                alpha = j.at("alpha");
                table = sf.generate( [&](double q) { return 1 - std::erfc(alpha*rc)*q + q*q; }, 0, 1 );
                calcDielectric = [&](double M2V) { return 1 + 3*M2V; };
                selfenergy_prefactor = erf(alpha*rc);
            }

            void sfFanourgakis(const json &j) {
                table = sf.generate( [&](double q) { return 1 - 1.75*q + 5.25*pow(q,5) - 7*pow(q,6) + 2.5*pow(q,7); }, 0, 1 );
                calcDielectric = [&](double M2V) { return 1 + 3*M2V; };
                selfenergy_prefactor = 0.875;
            }

            void sfFennel(const json &j) {
                alpha = j.at("alpha");
                table = sf.generate( [&](double q) { return (erfc(alpha*rc*q) - std::erfc(alpha*rc)*q + (q-1.0)*q*(std::erfc(alpha*rc)
                                + 2 * alpha * rc / std::sqrt(pc::pi) * std::exp(-alpha*alpha*rc*rc))); }, 0, 1 );
                calcDielectric = [&](double M2V) { double T = erf(alpha*rc) - (2 / (3 * sqrt(pc::pi)))
                    * exp(-alpha*alpha*rc*rc) * (alpha*alpha*rc*rc * alpha*alpha*rc*rc + 2.0 * alpha*alpha*rc*rc + 3.0);
                    return (((T + 2.0) * M2V + 1.0)/ ((T - 1.0) * M2V + 1.0)); };
                selfenergy_prefactor = ( erfc(alpha*rc)/2.0 + alpha*rc/sqrt(pc::pi) );
            }

            void sfEwald(const json &j) {
                alpha = j.at("alpha");
                table = sf.generate( [&](double q) { return std::erfc(alpha*rc*q); }, 0, 1 );
                calcDielectric = [&](double M2V) {
                    double T = std::erf(alpha*rc) - (2 / (3 * sqrt(pc::pi)))
                        * std::exp(-alpha*alpha*rc*rc) * ( 2*alpha*alpha*rc*rc + 3);
                    return ((T + 2.0) * M2V + 1)/ ((T - 1) * M2V + 1);
                };
                selfenergy_prefactor = ( erfc(alpha*rc) + alpha*rc/sqrt(pc::pi)*(1.0 + std::exp(-alpha*alpha*rc2)) );
            }

            void sfWolf(const json &j) {
                alpha = j.at("alpha");
                table = sf.generate( [&](double q) { return (erfc(alpha*rc*q) - erfc(alpha*rc)*q); }, 0, 1 );
                calcDielectric = [&](double M2V) { double T = erf(alpha*rc) - (2 / (3 * sqrt(pc::pi))) * exp(-alpha*alpha*rc*rc)
                    * ( 2.0 * alpha*alpha*rc*rc + 3.0);
                    return (((T + 2.0) * M2V + 1.0)/ ((T - 1.0) * M2V + 1.0));};
                selfenergy_prefactor = ( erfc(alpha*rc) + alpha*rc/sqrt(pc::pi)*(1.0 + exp(-alpha*alpha*rc2)) );
            }

            void sfPlain(const json &j, double val=1) {
                table = sf.generate( [&](double q) { return val; }, 0, 1 );
                calcDielectric = [&](double M2V) { return (2.0*M2V + 1.0)/(1.0 - M2V); };
                selfenergy_prefactor = 0.0;
            }

            public:
            CoulombGalore(const std::string &name="coulomb") { PairPotentialBase::name=name; }

            void from_json(const json &j) override {
                try {
                    type = j.at("type");
                    rc = j.at("cutoff");
                    rc2 = rc*rc;
                    rc1i = 1/rc;
                    epsr = j.at("epsr");
                    lB = pc::lB( epsr );
                    depsdt = j.value("depsdt", -0.368*pc::temperature/epsr);
                    sf.setTolerance(
                            j.value("utol",1e-5),j.value("ftol",1e-2) );

                    if (type=="reactionfield") sfReactionField(j);
                    if (type=="fanourgakis") sfFanourgakis(j);
                    if (type=="qpotential") sfQpotential(j);
                    if (type=="yonezawa") sfYonezawa(j);
                    if (type=="yukawa") sfYukawa(j);
                    if (type=="fennel") sfFennel(j);
                    if (type=="plain") sfPlain(j,1);
                    if (type=="ewald") sfEwald(j);
                    if (type=="none") sfPlain(j,0);
                    if (type=="wolf") sfWolf(j);
                    if ( table.empty() )
                        throw std::runtime_error(name + ": unknown coulomb type '" + type + "'" );
                }

                catch ( std::exception &e ) {
                    std::cerr << "CoulombGalore error: " << e.what();
                    throw;
                }
            }

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

            double dielectric_constant(double M2V) {
                return calcDielectric( M2V );
            } 

            void to_json(json &j) const override {
                using namespace u8;
                j["epsr"] = epsr;
                j["T"+partial+epsilon_m + "/" + partial + "T"] = depsdt;
                j["lB"] = lB;
                j["cutoff"] = rc;
                j["type"] = type;
                if (type=="yukawa") {
                    j["debyelength"] = 1.0/kappa;
                    j["ionic strength"] = I;
                }
                if (type=="qpotential")
                    j["order"] = order;
                if (type=="yonezawa" || type=="fennel" || type=="wolf" || type=="ewald")
                    j["alpha"] = alpha;
                if (type=="reactionfield") {
                    if(epsrf > 1e10)
                        j[epsilon_m+"_rf"] = 2e10;
                    else
                        j[epsilon_m+"_rf"] = epsrf;
                }
                _roundjson(j, 5);
            }
        };

        /**
         * @brief Arbitrary potentials for specific atom types
         *
         * This maintains a species x species matrix with function pointers (`std::function`)
         * that wraps pair potentials. Flexibility over performance.
         *
         * @todo `to_json` should retrive info from potentials instead of merely passing input
         * @warning Each atom pair will be assigned an instance of a pair-potential. This may be
         *          problematic is these have large memory requirements (Lennard-Jones for instance
         *          keep a pair-matrix with sigma/epsilon values. Make these static?)
         */
        template<class T /** particle type */>
            class FunctorPotential : public PairPotentialBase {
                typedef std::function<double(const T&, const T&, const Point&)> uFunc;
                PairMatrix<uFunc,true> umatrix; // matrix with potential for each atom pair
                json _j; // storage for input json

                uFunc combineFunc(const json &j) const {
                    uFunc u = [](const T&a, const T&b, const Point &r){return 0.0;};
                    if (j.is_array()) {
                        for (auto &i : j) // loop over all defined potentials in array
                            if (i.is_object() && (i.size()==1))
                                for (auto it=i.begin(); it!=i.end(); ++it) {
                                    uFunc _u = nullptr;
                                    if (it.key()=="coulomb") _u = CoulombGalore() = i;
                                    if (it.key()=="cos2") _u = CosAttract() = i;
                                    if (it.key()=="polar") _u = Polarizability<T>() = i;
                                    if (it.key()=="hardsphere") _u = HardSphere<T>() = i;
                                    if (it.key()=="lennardjones") _u = LennardJones<T>() = i;
                                    if (it.key()=="repulsionr3") _u = RepulsionR3() = i;
                                    if (it.key()=="wca") _u = WeeksChandlerAndersen<T>() = i;
                                    if (it.key()=="pm") _u = Coulomb() + HardSphere<T>() = it.value();
                                    if (it.key()=="pmwca") _u = Coulomb() + WeeksChandlerAndersen<T>() = it.value();
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
                    return umatrix(a.id, b.id)(a, b, r);
                }

                void to_json(json &j) const override { j = _j; }

                void from_json(const json &j) override {
                    _j = j;
                    umatrix = decltype(umatrix)( atoms<T>.size(), combineFunc(j.at("default")) );
                    for (auto it=j.begin(); it!=j.end(); ++it) {
                        auto atompair = words2vec<std::string>(it.key()); // is this for a pair of atoms?
                        if (atompair.size()==2) {
                            auto ids = names2ids(atoms<T>, atompair);
                            umatrix.set(ids[0], ids[1], combineFunc(it.value()));
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

                atoms<T> = j["atomlist"].get<decltype(atoms<T>)>();

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

            T a = atoms<T>[0].p;
            T b = atoms<T>[1].p;
            T c = atoms<T>[2].p;
            Point r={2,0,0};
            CHECK( u(a,a,r) == Approx( coulomb(a,a,r) ) );
            CHECK( u(b,b,r) == Approx( coulomb(b,b,r) ) );
            CHECK( u(a,b,r) == Approx( coulomb(a,b,r) + wca(a,b,r) ) );
            CHECK( u(c,c,r*1.01) == 0 );
            CHECK( u(c,c,r*0.99) == pc::infty );
        }
#endif

        struct BondData {
            enum Variant {harmonic, fene, yukawa, dihedral, none};
            Variant type=none;
            std::vector<int> index;
            std::vector<double> k;

            void shift( int offset ) {
                for ( auto &i : index )
                    i += offset;
            } // shift index

            template<class Tpvec>
                double energy(const Tpvec &p, Geometry::DistanceFunction dist) const {
                    double d, wca, x, zz;
                    switch (type) {
                        case BondData::harmonic:
                            d = k[1] - dist(p[index[0]].pos, p[index[1]].pos).norm();
                            return 0.5 * k[0]*d*d;
                        case BondData::fene:
                            d = dist( p[index[0]].pos, p[index[1]].pos ).squaredNorm();
                            wca = 0;
                            x = k[3];
                            if (d<=x*1.2599210498948732) {
                                x = x/d;
                                x = x*x*x;
                                wca = k[2]*(x*x - x + 0.25);
                            }
                            return (d>k[1]) ? pc::infty : -0.5*k[0]*k[1]*std::log(1-d/k[1]) + wca;
                        case BondData::yukawa:
                            d = dist( p[index[0]].pos, p[index[1]].pos ).norm();
                            zz = p[index[0]].charge*p[index[1]].charge;
                            if (zz>1e-5) 
                                return zz*k[0]*exp(-d/k[1])/d;
                            else
                                return 0;
                        default: break;
                    }
                    assert(!"not implemented");
                    return 0;
                } 
        }; //!< Harmonic and angular potentials for bonded interactions

        void to_json(json &j, const BondData &b) {
            switch(b.type) {
                case BondData::harmonic:
                    j["harmonic"] = {
                        { "index", b.index },
                        { "k", b.k[0] / 1.0_kJmol * std::pow(1.0_angstrom, 2) },
                        { "req", b.k[1] / 1.0_angstrom } };
                    break;
                case BondData::fene:
                    j["fene"] = {
                        { "index", b.index },
                        { "k", b.k[0] / (1.0_kJmol / std::pow(1.0_angstrom, 2)) },
                        { "rmax", std::sqrt(b.k[1]) / 1.0_angstrom },
                        { "eps", b.k[2] / 1.0_kJmol },
                        { "sigma", std::sqrt(b.k[3]) / 1.0_angstrom } };
                    break;
                case BondData::yukawa:
                    j["yukawa"] = {
                        { "index", b.index },
                        { "bjerrumlength", b.k[0] / 1.0_angstrom },
                        { "debyelength", b.k[1] / 1.0_angstrom } };
                    break;

                default: break;
            }
        }

        void from_json(const json &j, BondData &b) {
            if (j.is_object())
                if (j.size()==1) {
                    std::string t = j.begin().key();
                    auto& val = j.begin().value();
                    if (t=="harmonic") {
                        b.type = BondData::harmonic;
                        b.index = val.at("index").get<decltype(b.index)>();
                        if (b.index.size()!=2)
                            throw std::runtime_error("harmonic bond requires exactly two indeces");
                        b.k.resize(2);
                        b.k[0] = val.at("k").get<double>() * 1.0_kJmol / std::pow(1.0_angstrom, 2); // k
                        b.k[1] = val.at("req").get<double>() * 1.0_angstrom; // req
                        return;
                    }
                    if (t=="fene") {
                        b.type = BondData::fene;
                        b.index = val.at("index").get<decltype(b.index)>();
                        if (b.index.size()!=2)
                            throw std::runtime_error("FENE bond requires exactly two indeces");
                        b.k.resize(4);
                        b.k[0] = val.at("k").get<double>() * 1.0_kJmol / std::pow(1.0_angstrom, 2); // k
                        b.k[1] = std::pow( val.at("rmax").get<double>() * 1.0_angstrom, 2); // rmax^2
                        b.k[2] = val.at("eps").get<double>() * 1.0_kJmol; // epsilon wca
                        b.k[3] = std::pow( val.at("sigma").get<double>() * 1.0_angstrom, 2); // sigma^2 wca
                        return;
                    }
                    if (t=="yukawa") {
                        b.type = BondData::yukawa;
                        b.index = val.at("index").get<decltype(b.index)>();
                        if (b.index.size()!=2)
                            throw std::runtime_error("yukawa requires exactly two indeces");
                        b.k.resize(2);
                        b.k[0] = pc::lB(val.at("epsr").get<double>()) * 1.0_angstrom; // bjerrum length
                        b.k[1] = val.at("debyelength").get<double>() * 1.0_angstrom; // debye length
                        return;
                    }

                    if (t=="dihedral") {
                        b.type = BondData::dihedral;
                        assert(!"to be implemented");
                        return;
                    }
                    if (b.type==BondData::none)
                        throw std::runtime_error("unknown bondtype: " + t);
                    else return;
                }
            throw std::runtime_error("error parsing json to bond");
        }

        inline auto filterBonds(const std::vector<BondData> &bonds, BondData::Variant bondtype) {
            std::vector<std::reference_wrapper<const BondData>> filt;
            filt.reserve(bonds.size());
            std::copy_if(bonds.begin(), bonds.end(), std::back_inserter(filt),
                    [bondtype](const auto &d){return d.type==bondtype;} );
            return filt;
        } //!< Filter bond container for matching bond type and return _reference_ to original

#ifdef DOCTEST_LIBRARY_INCLUDED
        TEST_CASE("[Faunus] BondData")
        {
            BondData b;

            // test harmonic
            json j = R"( {"harmonic" : { "index":[2,3], "k":0.5, "req":2.1} } )"_json;
            b = j;
            CHECK( j == json(b) );
            CHECK_THROWS( b = R"( {"harmonic" : { "index":[2], "k":0.5, "req":2.1} } )"_json );
            CHECK_THROWS( b = R"( {"harmonic" : { "index":[2,3], "req":2.1} } )"_json );
            CHECK_THROWS( b = R"( {"harmonic" : { "index":[2,3], "k":2.1} } )"_json );

            // test fene
            j = R"( {"fene" : { "index":[2,3], "k":1, "rmax":2.1, "eps":2.48, "sigma":2 } } )"_json;
            b = j;
            CHECK( j == json(b) );
            CHECK_THROWS( b = R"( {"fene" : { "index":[2,3,4], "k":1, "rmax":2.1} } )"_json );
            CHECK_THROWS( b = R"( {"fene" : { "index":[2,3], "rmax":2.1} } )"_json );
            CHECK_THROWS( b = R"( {"fene" : { "index":[2,3], "k":1} } )"_json );
            CHECK_THROWS( b = R"( {"fene" : { "index":[2,3], "k":1, "rmax":2.1, "eps":2.48} } )"_json );
            CHECK_THROWS( b = R"( {"fene" : { "index":[2,3], "k":1, "rmax":2.1, "sigma":2} } )"_json );

            CHECK_THROWS( b = R"( {"unknown" : { "index":[2,3], "k":2.1, "req":1.0} } )"_json );
            j = json::object();
            CHECK_THROWS( b = j );

            BondData b2 = R"( {"harmonic" : { "index":[2,3], "k":0.5, "req":2.1} } )"_json;
            std::vector<BondData> bonds = {b,b2};
            auto filt = filterBonds(bonds, BondData::harmonic);
            CHECK( filt.size() == 1 ); 
            CHECK( filt[0].get().type == BondData::harmonic);
            CHECK( &filt[0].get() == &bonds[1] ); // filt should contain references to bonds
        }
#endif

    }//end of namespace Potential
}//end of namespace Faunus
