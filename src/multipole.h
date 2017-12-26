#pragma once
#include "tabulate.h"
#include "potentials.h"

namespace Faunus {
    /**
     * @brief Returns ion-dipole interaction.
     * @param QBxMuA Product of ion B:s charge and dipole A:s scalar
     * @param muA Unit dipole moment vector of particel A
     * @param QAxMuB Product of ion A:s charge and dipole B:s scalar
     * @param muB Unit dipole moment vector of particel B
     * @param r Direction \f$ r_A - r_B \f$
     */
    template<class Tvec>
        double q2mu(double QBxMuA, const Tvec &muA, double QAxMuB, const Tvec &muB, const Tvec &r) {
            double r2i = 1/r.squaredNorm();  // B = sol(dip), A = ch(charge)
            double r1i = sqrt(r2i);
            double r3i = r1i*r2i;
            double W1 = QBxMuA*muA.dot(r)*r3i;
            double W2 = QAxMuB*muB.dot(-r)*r3i;
            return (W1+W2);
        }

    /**
     * @brief Returns dipole-dipole interaction
     *
     * @param muA Unit dipole moment vector of particle A
     * @param muB Unit dipole moment vector of particle B
     * @param muAxmuB Product of dipole scalars
     * @param r Direction \f$ r_A - r_B \f$
     */
    template<class Tvec>
        double mu2mu(const Tvec &muA, const Tvec &muB, double muAxmuB, const Tvec &r, double a=1.0, double b=0.0) {
#ifdef FAU_APPROXMATH
            double r1i = invsqrtQuake( r.squaredNorm() );
            double r2i = r1i*r1i;
#else
            double r2i = 1/r.squaredNorm();
            double r1i = sqrt(r2i);
#endif
            double dot = muA.dot(muB);
            double T = (3*muA.dot(r)*muB.dot(r)*r2i - dot)*a + dot*b;
            return -muAxmuB*T*r1i*r2i;
        }

    /**
     * @brief Returns ion-quadrupole interaction
     * @param qA Charge of particle A
     * @param quadB Quadrupole moment of particle B
     * @param qB Charge of particle B
     * @param quadA Quadrupole moment of particle A
     * @param r Direction \f$ r_A - r_B \f$
     */
    template<class Tvec, class Tmat>
        double q2quad(double qA, const Tmat &quadB,double qB, const Tmat &quadA, const Tvec &r) {
            double r2i = 1/r.squaredNorm();
            double r1i = sqrt(r2i);
            double r3i = r1i*r2i;
            double r5i = r3i*r2i;
            double WAB = r.transpose()*quadB*r;
            WAB = 3*WAB*r5i - quadB.trace()*r3i;
            double WBA = r.transpose()*quadA*r;
            WBA = 3*WBA*r5i - quadA.trace()*r3i;
            return (qA*WAB + qB*WBA);
        }

    namespace Potential {

        /**
         * @brief Help-function for the q-potential in class CoulombGalore
         */
        inline double qPochhammerSymbol(double q, int k=1, int P=300) {
            //P = 300 gives an error of about 10^-17 for k < 4

            double value = 1.0;
            double temp = pow(q,k);
            for(int i = 0; i < P; i++) {
                value *= (1.0 - temp);
                temp *= q;
            }
            return value;
        }

        /**
         * @brief Coulomb type potentials with spherical cutoff
         *
         * Beyond a spherical cutoff, \f$R_c\f$, the potential is cut and if
         * below, \f$ u(r) = \frac{\lambda_B z_i z_j }{ r }\mathcal{S}(q) \f$ with \f$q=r/R_c\f$
         * is returned with the following splitting functions, \f$\mathcal{S}\f$, that
         * will be splined during construction and thus evaluate at similar speeds,
         *
         *  Type            | \f$\mathcal{S}(q=r/R_c)\f$               | Additional keywords  | Reference / Comment
         *  --------------- | ---------------------------------------- | -------------------- | ----------------------
         *  `plain`         | \f$ 1 \f$                                | none                 | http://doi.org/ctnnsj
         *  `none`          | \f$ 0 \f$                                | none                 | For convenience, only.
         *  `wolf`          | \f$ \text{erfc}(\alpha R_c q)-erfc(\alpha R_c)q \f$ | `alpha`          | http://doi.org/cfcxdk
         *  `fennel`        | \f$ erfc(\alpha R_c q)-erfc(\alpha R_c)q + ( q -1 ) q \left( erfc(\alpha R_c) + \frac{2\alpha R_c}{\sqrt{\pi}} e^{-\alpha^2 R_c^2} \right) \f$ | `alpha`| http://doi.org/bqgmv2
         *  `yonezawa`      | \f$ 1 + erfc(\alpha R_c)q + q^2 \f$      | `alpha`              | http://dx.doi.org/10/j97
         *  `fanourgakis`   | \f$ 1-\frac{7}{4}q+\frac{21}{4}q^5-7q^6+\frac{5}{2}q^7\f$| none | http://doi.org/f639q5
         *  `qpotential`     | \f$ \prod_{n=1}^{order}(1-q^n) \f$       | `order=300`          | ISBN [9789174224405](http://goo.gl/hynRTS) (Paper V)
         *  `reactionfield` | \f$ 1 + \frac{\varepsilon_{RF}-\varepsilon_{r}}{2\varepsilon_{RF}+\varepsilon_{r}} q^3  - 3\frac{\varepsilon_{RF}}{2\varepsilon_{RF}+\varepsilon_{r}}q \f$      | `epsrf`     | http://doi.org/dbs99w
         *  `yukawa`        | \f$ e^{-\kappa R_c q}-e^{-\kappa R_c}\f$  | `debyelength`        | ISBN 0486652424
         *
         *  The following keywords are required for all types:
         *
         *  Keyword       |  Description
         *  ------------- |  -------------------------------------------
         *  `coulombtype` |  Type of splitting function as defined above
         *  `cutoff`      |  Spherical cutoff in angstroms
         *  `epsr`        |  Relative dielectric constant of the medium
         *
         *  More info:
         * 
         *  - On the dielectric constant, http://dx.doi.org/10.1080/00268978300102721
         *  - Generalized reaction field using ionic strength, http://dx.doi.org/10.1063/1.469273
         */
        class CoulombGalore : public PairPotentialBase {
            private:
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
                    table = sf.generate( [&](double q) { return 1 - erfc(alpha*rc)*q + q*q; }, 0, 1 );
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
                    table = sf.generate( [&](double q) { return (erfc(alpha*rc*q) - erfc(alpha*rc)*q + (q-1.0)*q*(erfc(alpha*rc)
                                    + 2 * alpha * rc / sqrt(pc::pi) * exp(-alpha*alpha*rc*rc))); }, 0, 1 );
                    calcDielectric = [&](double M2V) { double T = erf(alpha*rc) - (2 / (3 * sqrt(pc::pi)))
                        * exp(-alpha*alpha*rc*rc) * (alpha*alpha*rc*rc * alpha*alpha*rc*rc + 2.0 * alpha*alpha*rc*rc + 3.0);
                        return (((T + 2.0) * M2V + 1.0)/ ((T - 1.0) * M2V + 1.0)); };
                    selfenergy_prefactor = ( erfc(alpha*rc)/2.0 + alpha*rc/sqrt(pc::pi) );
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
                CoulombGalore() { name="coulomb"; }
                
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
                                j.value("utol",1e-9),j.value("ftol",1e-2) );

                        if (type=="reactionfield") sfReactionField(j);
                        if (type=="fanourgakis") sfFanourgakis(j);
                        if (type=="qpotential") sfQpotential(j);
                        if (type=="yonezawa") sfYonezawa(j);
                        if (type=="yukawa") sfYukawa(j);
                        if (type=="fennel") sfFennel(j);
                        if (type=="plain") sfPlain(j,1);
                        if (type=="none") sfPlain(j,0);
                        if (type=="wolf") sfWolf(j);
                        if ( table.empty() )
                            throw std::runtime_error(name + ": unknown coulomb type '" + type + "'" );
                    }

                    catch ( std::exception &e )
                    {
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
                        j["debye length"] = 1.0/kappa;
                        j["ionic strength"] = I;
                    }
                    if (type=="qpotential")
                        j["order"] = order;
                    if (type=="yonezawa" || type=="fennel" || type=="wolf")
                        j["alpha"] = alpha;
                    if (type=="reactionfield") {
                        if(epsrf > 1e10)
                            j[epsilon_m+"_rf"] = pc::infty;
                        else
                            j[epsilon_m+"_rf"] = epsrf;
                    }
                    _roundjson(j, 5);
                }
        };

    } // potential namespace
} // faunus namespace

