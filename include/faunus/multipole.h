#ifndef FAUNUS_MULTIPOLE_H
#define FAUNUS_MULTIPOLE_H
#define DIPOLEPARTICLE

#include <faunus/common.h>
#include <faunus/auxiliary.h>
#include <faunus/species.h>
#include <faunus/tabulate.h>

namespace Faunus {
    namespace json {
        /**
         * @brief Loads a JSON file, reads atom pair properties and returns a vector map
         *
         * Example:
         * ~~~~
         * auto map = atomPairMap("input.json", "pairproperties", "nemorep");
         * for (auto &m : map)
         *   cout << m.second.transpose() << endl; // -> 12 23 0.2 -2 3 4 5 ...
         * ~~~~
         * where the input json file could look like this:
         * ~~~~
         * {
         *   "pairproperties" : {
         *      "OW OW"  : { "nemorep":"12. 23. 0.2  -2   3   4   5" },
         *      "HW HW"  : { "nemorep":"-2. 23. 0.2   2  99   4  -5 " },
         *      "HW OW"  : { "nemorep":"112. 23. 0.2 129 391 238  23" }
         *   }
         * }
         * ~~~~
         */
        inline std::map<opair<int>,Eigen::VectorXd>
            atomPairMap(const string &file, const string &section, const string &key) {
                assert(!section.empty() && !key.empty());
                typedef Eigen::VectorXd Tvec;
                typedef opair<int> Tpair;
                std::map<Tpair,Tvec> map;
                string atom1, atom2;
                /*
                   auto j=json::open(file);
                   for (auto &a : json::object(section, j)) {
                   std::istringstream is(a.first);
                   is >> atom1 >> atom2;
                   Tpair pair( atom[atom1].id, atom[atom2].id );
                   string str = json::value<string>(a.second, key,"");
                   std::istringstream is2(str), tmp(str);
                   int size = std::distance(std::istream_iterator<string>(tmp), std::istream_iterator<string>());
                   Tvec v(size);
                   for (int i=0; i<size; i++)
                   is2 >> v[i];
                   map[pair] = v;
                   }*/
                return map;
            }
    }//namespace

    /**
     * @brief Approximation of erfc-function
     * @param x Value for which erfc should be calculated 
     * @details Reference for this approximation is found in Abramowitz and Stegun, 
     *          Handbook of mathematical functions, eq. 7.1.26
     *
     * @f[
     *     erf(x) = 1 - (a_1t + a_2t^2 + a_3t^3 + a_4t^4 + a_5t^5)e^{-x^2} + \epsilon(x)
     * @f]
     * @f[
     *     t = \frac{1}{1 + px}
     * @f]
     * @f[
     *     |\epsilon(x)| \le 1.5\times 10^{-7}
     * @f]
     * 
     * @warning Needs modification if x < 0
     */
    inline double erfc_x(double x) {
        // |error| <= 1.5*10^-7
        if(x < 0.0)
            return ( 2.0 - erfc_x(-x) );

        double p = 0.3275911;
        double t = 1.0/(1.0+p*x);
        double x2 = x*x;
        double a1 = 0.254829592;
        double a2 = -0.284496736;
        double a3 = 1.421413741;
        double a4 = -1.453152027;
        double a5 = 1.061405429;
        double tp = t*(a1+t*(a2+t*(a3+t*(a4+t*a5))));
        return tp*exp(-x2);
    }

    /**
     * @brief See erfc_x-function. 1 - erfc_x
     * @param x Value for which erf should be calculated 
     */
    inline double erf_x(double x) {
        return (1 - erfc_x(x));
    }

    /**
     * @brief Returns NemoType1-interaction (Exponential Repulsion)                         Needs to be checked!
     * @param vec Vector with parameters. Form: \f$ [a_{ab} b_{ab} c_{ab} d_{ab}] \f$  
     * @param r Vector between particles
     * @param expmax Maximum exponential coefficient (optional)
     */
    template<class Tvec>
        double nemo1(Eigen::VectorXd &vec, const Tvec &r, double expmax=80.0) {
            double asw = 1.2;
            double nsw = 4;
            double r1i = 1/r.norm();
            double r2i = r1i*r1i;
            double r6i = r2i*r2i*r2i;
            double ss = 1 - exp(-std::min(expmax,pow((1/(asw*r1i)),nsw)));

            double uexp  = vec[0]*exp(-std::min(expmax,vec[1]/r1i));
            double ur20  = vec[2]*r6i*r6i*r6i*r2i;
            double udis  = vec[3]*ss*r6i;
            return (uexp  + ur20 + udis);
        }

    /**
     * @brief Returns NemoType2-interaction (r-7 Repulsion)                         Needs to be checked!
     * @param vec Vector with parameters. Form: \f$ [a_{ab} b_{ab} c_{ab} d_{ab}] \f$  
     * @param r Vector between particles
     * @param expmax Maximum exponential coefficient (optional)
     */
    template<class Tvec>
        double nemo2(Eigen::VectorXd &vec, const Tvec &r, double expmax=80.0) {
            double asw = 1.2;
            double nsw = 4;
            double r1i = 1/r.norm();
            double r2i = r1i*r1i;
            double r6i = r2i*r2i*r2i;
            double ss = 1 - exp(-std::min(expmax,pow((1/(asw*r1i)),nsw)));

            double uexp  = vec[0]*r1i*r6i;
            double udis  = vec[3]*ss*r6i;
            return (uexp + udis);
        }

    /**
     * @brief Returns NemoType3-interaction (Modified Interactions)                         Needs to be checked!
     * @param vec Vector with parameters. Form: \f$ [a_{ab} b_{ab} c_{ab} d_{ab} n_{ab}] \f$  
     * @param r Vector between particles
     */
    template<class Tvec>
        double nemo3(Eigen::VectorXd &vec, const Tvec &r) {
            double r1i = 1/r.norm();
            double r2i = r1i*r1i;
            double r6i = r2i*r2i*r2i;

            double uexp  = vec[3]*pow(r1i,vec[4]);          // vec[4] = nab   <------------------------
            double udis1  = -vec[2]*r6i;
            double udis2  = vec[0]*exp(-vec[1]/r1i);
            return (uexp  + udis1 + udis2);
        }

    /**
     * @brief Returns NemoType4-interaction (Damping Exponential)
     * @param vec Vector with parameters. Form: \f$ [a_{ab} b_{ab} c_{ab} d_{ab} e_{ab} f_{ab} n_{ab}] \f$  
     * @param r Vector between particles
     * @param expmax Maximum exponential coefficient (optional)
     */
    template<class Tvec>
        double nemo4(Eigen::VectorXd &vec, const Tvec &r, double expmax=80.0) {
            double r1i = 1/r.norm();
            double r2i = r1i*r1i;
            double r6i = r2i*r2i*r2i;

            double uexp1  = vec[4]*exp(-std::min(expmax,vec[5]/r1i));
            double uexp2  = 0;
            if(vec[6] != 0) uexp2 = vec[3]*pow(r1i,vec[6]);               // vec[6] = nab   <------------------------
            double udis1  =-vec[2]*r6i;
            double udis2  = vec[0]*exp(-std::min(expmax,vec[1]/r1i));
            return (uexp1 + uexp2  + udis1 + udis2);
        }

    /**
     * @brief Returns NemoType5-interaction (Full Damping)                         Needs to be checked!
     * @param vec Vector with parameters. Form: \f$ [a_{ab} b_{ab} c_{ab} d_{ab} e_{ab} f_{ab} n_{ab}] \f$  
     * @param r Vector between particles
     * @param expmax Maximum exponential coefficient (optional)
     */
    template<class Tvec>
        double nemo5(Eigen::VectorXd &vec, const Tvec &r, double expmax=80.0) {
            double r1i = 1/r.norm();
            double r2i = r1i*r1i;
            double r6i = r2i*r2i*r2i;
            double bri = r1i/vec[1];
            double ud1 = 6*bri;
            double ud2 = 5*bri*ud1;
            double ud3 = 4*bri*ud2;
            double ud4 = 3*bri*ud3;
            double ud5 = 2*bri*ud4;
            double ud6 = bri*ud5;

            double uexp1  = vec[4]*exp(-std::min(expmax,vec[5]/r1i));
            double uexp2  = 0;
            if(vec[6] != 0) uexp2 = vec[3]*pow(r1i,vec[6]);       // vec[6] = nab   <------------------------
            double udis1  =-vec[2]*r6i;
            double udd = 1 + ud1 + ud2 + ud3 + ud4 + ud5 + ud6;
            double udis2  = vec[0]*exp(-std::min(expmax,1/bri));
            return (uexp1 + uexp2  + udis1 + udd*udis2);
        }

    /**
     * @brief Returns NemoType6-interaction (Full Damping chtr)                         Needs to be checked!
     * @param vec Vector with parameters. Form: \f$ [a_{ab} b_{ab} c_{ab} d_{ab} e_{ab} f_{ab} n_{ab}] \f$  
     * @param r Vector between particles
     * @param expmax Maximum exponential coefficient (optional)
     */
    template<class Tvec>
        double nemo6(Eigen::VectorXd &vec, const Tvec &r, double expmax=80.0) {
            double r1i = 1/r.norm();
            double r2i = r1i*r1i;
            double r6i = r2i*r2i*r2i;
            double bri = r1i/vec[1];
            double ud1 = 6*bri;
            double ud2 = 5*bri*ud1;
            double ud3 = 4*bri*ud2;
            double ud4 = 3*bri*ud3;
            double ud5 = 2*bri*ud4;
            double ud6 = bri*ud5;

            double uexp1  = vec[4]*exp(-std::min(expmax,vec[5]/r1i));
            double uexp2  = 0;
            if(vec[6] != 0) uexp2 = vec[3]*pow(r1i,vec[6]);                     // vec[6] = nab   <------------------------
            double udis1  =-vec[2]*r6i;
            double udd = 1 + ud1 + ud2 + ud3 + ud4 + ud5 + ud6;
            double udis2  = vec[0]*exp(-std::min(expmax,1/bri));
            double uchtexp = -vec[8]*exp(-std::min(expmax,vec[7]/r1i));  // vec[7] = acht, vec[8] = kcht    <------------------------
            return (uexp1 + uexp2  + udis1 + udd*udis2 + uchtexp);
        }

    /**
     * @brief Returns NemoType7-interaction (Full Damping chtr gaussian)                         Needs to be checked!
     * @param vec Vector with parameters. Form: \f$ [a_{ab} b_{ab} c_{ab} d_{ab} e_{ab} f_{ab} n_{ab}  a_{cht} k_{cht}] \f$  
     * @param r Vector between particles
     * @param expmax Maximum exponential coefficient (optional)
     */
    template<class Tvec>
        double nemo7(Eigen::VectorXd &vec, const Tvec &r, double expmax=80.0) {
            double r1i = 1/r.norm();
            double r2i = r1i*r1i;
            double r6i = r2i*r2i*r2i;
            double bri = r1i/vec[1];
            double ud1 = 6*bri;
            double ud2 = 5*bri*ud1;
            double ud3 = 4*bri*ud2;
            double ud4 = 3*bri*ud3;
            double ud5 = 2*bri*ud4;
            double ud6 = bri*ud5;

            double uchtexp = -vec[8]*exp(-std::min(expmax,vec[7]*(pow((r.norm()-vec[3]),2))));  // vec[7] = acht, vec[8] = kcht    <------------------------
            double uexp  = vec[4]*exp(-std::min(expmax,vec[5]/r1i));
            double udis1  =-vec[2]*r6i;
            double udd = 1 + ud1 + ud2 + ud3 + ud4 + ud5 + ud6;
            double udis2  = vec[0]*exp(-std::min(expmax,1/bri));
            return (uexp + udis1 + udd*udis2 + uchtexp);
        }

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

        class NemoRepulsion : public PairPotentialBase {
            private:
                string _brief() { return "NemoRepulsion"; }
            protected:
                typedef Eigen::VectorXd Tvec;
                typedef opair<int> Tpair;
                std::map<Tpair,Tvec> pairMap;
                double expmax;
                double scaling;

            public:
                NemoRepulsion(InputMap &in) {
                    name="Nemo repulsion";
                    pc::setT ( in.get<double>("temperature", 298.15, "Absolute temperature (K)") );
                    expmax = in.get<double>("expmax", 80, "Maximum repulsion exponent");
                    scaling = 1000/(pc::Nav*pc::kT());  // Converts from kJ/mol to kT
                    pairMap = json::atomPairMap("water2.json","pairproperties","nemorep");
                }

                /**
                 * @brief NemoRepulsion
                 * @param a Dipole particle A
                 * @param b Dipole particle B
                 * @param r Direction \f$ r_A - r_B \f$  
                 */
                template<class Tparticle>
                    double operator()(const Tparticle &a, const Tparticle &b, const Point &r) const {
                        Tpair pair(a.id,b.id);
                        Tvec vec;
                        auto it = pairMap.find(pair);
                        if (it!=pairMap.end()) { 
                            vec = it->second;
                            return nemo4(vec, r,expmax)*scaling;
                        }
                        assert(!"No pair data defined");
                        return 0;
                    }

                string info(char w) { return _brief(); }
        };

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
         *  `wolf`          | \f$ erfc(\alpha R_c q)-erfc(\alpha R_c)q \f$ | `alpha`          | http://doi.org/cfcxdk
         *  `fennel`        | \f$ erfc(\alpha R_c q)-erfc(\alpha R_c)q + ( q -1 ) q \left( erfc(\alpha R_c) + \frac{2\alpha R_c}{\sqrt{\pi}} e^{-\alpha^2 R_c^2} \right) f$ | `alpha`| http://doi.org/bqgmv2
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
                string type;
		double selfenergy_prefactor;
                double lB, depsdt, rc, rc2, rc1i, epsr, epsrf, alpha, kappa, I;
                int order;

                void sfYukawa(const Tmjson &j) {
                    kappa = 1.0 / j.at("debyelength").get<double>();
                    I = kappa*kappa / ( 8.0*lB*pc::pi*pc::Nav/1e27 );
                    table = sf.generate( [&](double q) { return std::exp(-q*rc*kappa) - std::exp(-kappa*rc); } ); // q=r/Rc 
                    // we could also fill in some info string or JSON output...
                }

                void sfReactionField(const Tmjson &j) {
                    epsrf = j.at("eps_rf");
                    table = sf.generate( [&](double q) { return 1 + (( epsrf - epsr ) / ( 2 * epsrf + epsr ))*q*q*q - 3 * ( epsrf / ( 2 * epsrf + epsr ))*q ; } ); 
                    calcDielectric = [&](double M2V) {
                        if(epsrf > 1e10)
                            return 1 + 3*M2V;
                        if(fabs(epsrf-epsr) < 1e-6)
                            return 2.25*M2V + 0.25 + 0.75*sqrt(9*M2V*M2V + 2*M2V + 1);
                        if(fabs(epsrf-1.0) < 1e-6)
                            return ( 2*M2V + 1 ) / ( 1 - M2V );
                        return 0.5 * ( 2*epsrf - 1 + sqrt( -72*M2V*M2V*epsrf + 4*epsrf*epsrf + 4*epsrf + 1) ) / ( 3*M2V-1 ); // Needs to be checked!
                        //return (6*M2V*epsrf + 2*epsrf + 1.0)/(1.0 + 2*epsrf - 3*M2V); // Is OK when epsr=1.0
                    };
		    selfenergy_prefactor = 1.5*epsrf/(2.0*epsrf + epsr); // Correct?!, see Eq.14 in DOI: 10.1021/jp510612w
                    // we could also fill in some info string or JSON output...
                }

                void sfQpotential(const Tmjson &j)
                {
                    order = j.value("order",300);
                    table = sf.generate( [&](double q) { return qPochhammerSymbol( q, 1, order ); } );
                    calcDielectric = [&](double M2V) { return 1 + 3*M2V; };
		    selfenergy_prefactor = 0.5;
                }

                void sfYonezawa(const Tmjson &j)
                {
                    alpha = j.at("alpha");
                    table = sf.generate( [&](double q) { return 1 - erfc(alpha*rc)*q + q*q; } );
		    calcDielectric = [&](double M2V) { return 1 + 3*M2V; };
		    selfenergy_prefactor = erf(alpha*rc);
                }

                void sfFanourgakis(const Tmjson &j) {
                    table = sf.generate( [&](double q) { return 1 - 1.75*q + 5.25*pow(q,5) - 7*pow(q,6) + 2.5*pow(q,7); } );
                    calcDielectric = [&](double M2V) { return 1 + 3*M2V; };
		    selfenergy_prefactor = 0.875;
                }

                void sfFennel(const Tmjson &j) {
                    alpha = j.at("alpha");
                    table = sf.generate( [&](double q) { return (erfc(alpha*rc*q) - erfc(alpha*rc)*q + (q-1.0)*q*(erfc(alpha*rc) + 2 * alpha * rc / sqrt(pc::pi) * exp(-alpha*alpha*rc*rc))); } );
		    calcDielectric = [&](double M2V) { double T = erf(alpha*rc) - (2 / (3 * sqrt(pc::pi))) * exp(-alpha*alpha*rc*rc) * (alpha*alpha*rc*rc * alpha*alpha*rc*rc + 2.0 * alpha*alpha*rc*rc + 3.0);
						       return (((T + 2.0) * M2V + 1.0)/ ((T - 1.0) * M2V + 1.0)); };
		    selfenergy_prefactor = ( erfc(alpha*rc)/2.0 + alpha*rc/sqrt(pc::pi) );
                }

                void sfWolf(const Tmjson &j) {
                    alpha = j.at("alpha");
                    table = sf.generate( [&](double q) { return (erfc(alpha*rc*q) - erfc(alpha*rc)*q); } );
		    calcDielectric = [&](double M2V) { double T = erf(alpha*rc) - (2 / (3 * sqrt(pc::pi))) * exp(-alpha*alpha*rc*rc) * ( 2.0 * alpha*alpha*rc*rc + 3.0);
						       return (((T + 2.0) * M2V + 1.0)/ ((T - 1.0) * M2V + 1.0));};
		    selfenergy_prefactor = ( erfc(alpha*rc) + alpha*rc/sqrt(pc::pi)*(1.0 + exp(-alpha*alpha*rc2)) );
                }

                void sfPlain(const Tmjson &j, double val=1) {
                    table = sf.generate( [&](double q) { return val; } );
		    calcDielectric = [&](double M2V) { return (2.0*M2V + 1.0)/(1.0 - M2V); };
		    selfenergy_prefactor = 0.0;
                }

            public:
                CoulombGalore(const Tmjson &j) {
                    try {
                        type = j.at("coulombtype");
                        name = "Coulomb-" + textio::toupper_first( type );
                        rc = j.at("cutoff");
                        rc2 = rc*rc;
                        rc1i = 1/rc;
                        epsr = j.at("epsr");
                        lB = pc::lB( epsr );
                        depsdt = j.value("depsdt", -0.368*pc::T()/epsr);

                        sf.setRange(0, 1);
                        sf.setTolerance(
                                j.value("tab_utol",1e-9),j.value("tab_ftol",1e-2) );

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
                            throw std::runtime_error("unknown coulomb type '" + type + "'" );
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
                            double r = sqrt(r2);
                            return lB * a.charge * b.charge / r * sf.eval( table, r*rc1i );
                        }
                        return 0;
                    }

                template<class Tparticle>
                    double operator()(const Tparticle &a, const Tparticle &b, const Point &r) const {
                        return operator()(a,b,r.squaredNorm());
                    }

                template<typename Tparticle>
                    Point force(const Tparticle &a, const Tparticle &b, double r2, const Point &p) {
                        if (r2 < rc2) {
                            double r = sqrt(r2);
                            return lB * a.charge * b.charge * ( -sf.eval( table, r*rc1i )/r2 + sf.evalDer( table, r*rc1i )/r )*p;
                        }
                        return Point(0,0,0);
                    }
                    
                /**
		 * @brief Self-energy of the potential
		 */
		template<class Tpvec>
		    double internal(const Tpvec &p, const Group &g) const { 
		    double Eq = 0;
		    for (auto i : g)
			Eq += p[i].charge * p[i].charge;
		    return -selfenergy_prefactor*Eq*lB/rc;
		}

                double dielectric_constant(double M2V) {
                    return calcDielectric( M2V );
                } 

                string info(char w) {
                    using namespace textio;
                    std::ostringstream o;
                    o << pad(SUB, w, "Temperature") << pc::T() << " K" << endl
                        << pad(SUB, w, "Dielectric constant") << epsr << endl
                        << pad(SUB, w + 6, "T" + partial + epsilon + "/" + epsilon + partial + "T") << depsdt << endl
                        << pad(SUB, w, "Bjerrum length") << lB << _angstrom << endl
                        << pad(SUB,w,"Cutoff") << rc << _angstrom << endl;
                    if (type=="yukawa") {
                        o << pad(SUB,w, "Debye length") << 1.0/kappa << endl;
                        o << pad(SUB,w, "Ionic strength") << I << endl;
                    }
                    if (type=="reactionfield") {
                        if(epsrf > 1e10)
                            o << pad(SUB,w+1, epsilon_m+"_RF") << infinity << endl;
                        else
                            o << pad(SUB,w+1, epsilon_m+"_RF") << epsrf << endl;
                    }
                    if (type=="qpotential")
                        o << pad(SUB,w, "order") << order << endl;
                    if (type=="yonezawa" || type=="fennel" || type=="wolf")
                        o << pad(SUB,w, "alpha") << alpha << endl;
                    return o.str();
                }
        };

        /**
         * @brief Help-function for `sfQpotential` using order 3.
         */
        inline void _DipoleDipoleQ2Help_3(double &a3, double &b3, double q) {
            double q2 = q*q;
            a3 = ( ( (-5.0 * q2 + 8.0 / 3.0 * q)*q + q)*q + 1.0 / 3.0)*q2 + 1.0;
            b3 = -2.0 / 3.0 * (  (15.0 * q2 - 10.0 * q - 5.0 )*q2 + 1.0) * q2;
        }

        /**
         * @brief Help-function for `sfQpotential` using order 4.
         */
        inline void _DipoleDipoleQ2Help_4(double &a4, double &b4, double q) {
            double q2 = q*q;
            double q3 = q2*q;
            a4 = ( ( (10.0 * q2 - 9.0 * q - 8.0 )*q3 + 10.0 )*q3 - 2.0 )*q - 1.0;
            b4 = ( ( ( 21.0 * q2 - 16.0 * q - 35.0 / 3.0 )*q3 + 16.0 / 3.0 )*q3 + 1.0 / 3.0 )*q2 + 1.0;
        }

        /**
         * @brief Help-function for `sfQpotential`.
         */
        inline double _DipoleDipoleQ2Help(double q, int l=0, int P=300, bool all=true) {
            if(q >= 1.0 - (1.0/2400.0))
                return 0.0;
            if(q <= (1.0/2400.0) && all)
                return 1.0;
            if(q <= (1.0/2400.0) && !all)
                return 0.0;

            double Q = 0.0;
            double T = 0.0;
            double qP = 1.0; // Will end as q-Pochhammer Symbol, (q^l;q)_P
            double fac = pow(q,l);
            for( int n = 1; n <= P; n++) {
                double temp2 = (n + l)/(1.0 - pow(q,(n+l)));
                fac *= q;
                Q -= temp2*fac;
                T -= temp2*temp2*fac;
                qP *= (1.0 - fac);
            }

            if(!all)
                return qP*(-Q + T + Q*Q)/3.0;
            return qP*(-Q + T + Q*Q - 3.0*Q + 3.0)/3.0;
        }
	
        /**
         * @brief Dipole-dipole type potentials with spherical cutoff
         *
         * Beyond a spherical cutoff, \f$R_c\f$, the potential is cut and if
         * below, the modified interaction tensor \f$ {\bf T}^{\rm Mod}_2 =  {\bf T}_2\mathcal{S}_1(q) + {\bf I}\mathcal{S}_2(q) \f$ is used where \f$ {\bf T}_2 \f$ 
	 * is the original dipole-dipole interaction-tensor, \f$q=r/R_c\f$, and \f$\mathcal{S}_1\f$ and \f$\mathcal{S}_2\f$ are splitting functions. 
	 * The interactions will thus be evaluated at similar speeds,
         *
         *  Type            | Additional keywords  | Reference / Comment
         *  --------------- | -------------------- | ----------------------
         *  `plain`         | none                 | http://doi.org/ctnnsj
         *  `none`          | none                 | For convenience, only.
         *  `wolf`          | `alpha`              | http://doi.org/cfcxdk
         *  `fennel`        | `alpha`              | http://doi.org/bqgmv2
         *  `fanourgakis`   |                      | http://doi.org/f639q5
         *  `qpotential`    | `order=300`          | 
	 *  `q2potential`   | `order=300`          | 
         *  `reactionfield` | `epsrf`              | http://doi.org/dbs99w
         *
         *  The following keywords are required for all types:
         *
         *  Keyword       |  Description
         *  ------------- |  -------------------------------------------
         *  `coulombtype` |  Type of splitting function as defined in reference above
         *  `cutoff`      |  Spherical cutoff in angstroms
         *  `epsr`        |  Relative dielectric constant of the medium
         *
         *  More info:
         * 
         *  - On the dielectric constant, http://dx.doi.org/10.1080/00268978300102721
         */
        class DipoleDipoleGalore : public PairPotentialBase {
            private:
                Tabulate::Andrea<double> ak, bk; // splitting function
                Tabulate::TabulatorBase<double>::data tableA, tableB; // data for splitting function
                std::function<double(double)> calcDielectric; // function for dielectric const. calc.
                string type;
		double selfenergy_prefactor;
                double lB, depsdt, rc, rc2, rc1i, epsr, epsrf, alpha;
		int order;

                void sfReactionField(double epsrf_in) {
                    epsrf = epsrf_in;
                    tableA = ak.generate( [&](double q) { return 1.0; } ); 
		    tableB = bk.generate( [&](double q) { return -(2*(epsrf-epsr)/(2*epsrf+epsr))/epsr*q*q*q; } ); 
                    calcDielectric = [&](double M2V) {
                        if(epsrf > 1e10)
                            return 1 + 3*M2V;
                        if(fabs(epsrf-epsr) < 1e-6)
                            return 2.25*M2V + 0.25 + 0.75*sqrt(9*M2V*M2V + 2*M2V + 1);
                        if(fabs(epsrf-1.0) < 1e-6)
                            return ( 2*M2V + 1 ) / ( 1 - M2V );
                        return 0.5 * ( 2*epsrf - 1 + sqrt( -72*M2V*M2V*epsrf + 4*epsrf*epsrf + 4*epsrf + 1) ) / ( 3*M2V-1 ); // Needs to be checked!
                        //return (6*M2V*epsrf + 2*epsrf + 1.0)/(1.0 + 2*epsrf - 3*M2V); // Is OK when epsr=1.0
			selfenergy_prefactor = 2.0*(epsr - epsrf)/(2.0*epsrf + epsr); // Preliminary, needs to be checked!
                    };
                }

                void sfQ2potential(int order_in)
                {
                    order = order_in;
                    tableA = ak.generate( [&](double q) { return qPochhammerSymbol(q,3,order); } );
                    tableB = bk.generate( [&](double q) { return 0.0; } );
		    calcDielectric = [&](double M2V) { return (2*M2V + 1.0)/(1.0 - M2V); };
		    selfenergy_prefactor = 0.5;
                }

                void sfQpotential(int order_in)
                {
                    order = order_in;
                    tableA = ak.generate( [&](double q) { return _DipoleDipoleQ2Help(q,0,order); } );
                    tableB = bk.generate( [&](double q) { return _DipoleDipoleQ2Help(q,0,order,false); } );
		    calcDielectric = [&](double M2V) { return 1 + 3*M2V; };
		    selfenergy_prefactor = 0.5;
                }

                void sfFanourgakis() {
                    tableA = ak.generate( [&](double q) { return ( 1.0 + 14.0*pow(q,5) - 35.0*pow(q,6) + 20.0*pow(q,7) ); } );
                    tableB = bk.generate( [&](double q) { return 35.0*pow(q,5)*( 1.0 - 2.0*q + q*q ); } );
		    calcDielectric = [&](double M2V) { return 1 + 3*M2V; };
		    selfenergy_prefactor = 0.0; // Seems so but is it really correct? Check!
                }

                void sfFennel(double alpha_in) {
                    alpha = alpha_in;
		    double ar = alpha*rc;
                    tableA = ak.generate( [&](double q) { return ( ( 2.0*ar*q*exp(-ar*ar*q*q)*(2.0*ar*ar*q*q/3.0 + 1.0)/sqrt(pc::pi) + erfc(ar*q)) - (erfc(ar) + 2.0*ar*exp(-ar*ar)*(2.0*ar*ar/3.0 + 1.0)/sqrt(pc::pi))*q*q*q + q*q*q*(q-1.0)*( 3.0*erfc(ar) + 2.0*ar*exp(-ar*ar)*(3.0 + 2.0*ar*ar + 4.0/3.0*ar*ar*ar*ar)/sqrt(pc::pi))); } );
		    tableB = bk.generate( [&](double q) { return ( 4.0/3.0*ar*ar*ar*q*q*q*exp(-ar*ar*q*q)/sqrt(pc::pi) - 4.0/3.0*ar*ar*ar*exp(-ar*ar)/sqrt(pc::pi)*q*q*q + q*q*q*(q - 1.0)*8.0/3.0*exp(-ar*ar)*pow(ar,5)/sqrt(pc::pi) ); } );
		    calcDielectric = [&](double M2V) { double T = erf(ar) - (2 / (3 * sqrt(pc::pi))) * exp(-ar*ar) * (ar*ar*ar*ar + 2.0 * ar*ar + 3.0);
						       return (((T + 2.0) * M2V + 1.0)/ ((T - 1.0) * M2V + 1.0)); };
		    selfenergy_prefactor = ( erfc(ar)/2.0 + ar/sqrt(pc::pi)*exp(-ar*ar) + (2.0/3.0)*ar*ar*ar/sqrt(pc::pi) );
                }

                void sfWolf(double alpha_in) {
                    alpha = alpha_in;
		    double ar = alpha*rc;
                    tableA = ak.generate( [&](double q) { return ( ( 2.0*ar*q*exp(-ar*ar*q*q)*(2.0*ar*ar*q*q/3.0 + 1.0)/sqrt(pc::pi) + erfc(ar*q)) - (erfc(ar) + 2.0*ar*exp(-ar*ar)*(2.0*ar*ar/3.0 + 1.0)/sqrt(pc::pi))*q*q*q ); } );
		    tableB = bk.generate( [&](double q) { return ( 4.0/3.0*ar*ar*ar*q*q*q*exp(-ar*ar*q*q)/sqrt(pc::pi) - 4.0/3.0*ar*ar*ar*exp(-ar*ar)/sqrt(pc::pi)*q*q*q ); } );
		    calcDielectric = [&](double M2V) { double T = erf(alpha*rc) - (2 / (3 * sqrt(pc::pi))) * exp(-alpha*alpha*rc*rc) * ( 2.0 * alpha*alpha*rc*rc + 3.0);
						       return (((T + 2.0) * M2V + 1.0)/ ((T - 1.0) * M2V + 1.0));};
		    selfenergy_prefactor = ( erfc(ar)/2.0 + ar/sqrt(pc::pi)*exp(-ar*ar) + (2.0/3.0)*ar*ar*ar/sqrt(pc::pi) );
                }

                void sfPlain(double val=1) {
                    tableA = ak.generate( [&](double q) { return val; } );
		    tableB = bk.generate( [&](double q) { return 0.0; } );
		    calcDielectric = [&](double M2V) { return (2.0*M2V + 1.0)/(1.0 - M2V); };
		    selfenergy_prefactor = 0.0;
                }

            public:
                DipoleDipoleGalore(const Tmjson &j) {
                    try {
                        type = j.at("coulombtype");
                        name = "DipoleDipole-" + textio::toupper_first( type );
                        rc = j.at("cutoff");
                        rc2 = rc*rc;
                        rc1i = 1/rc;
                        epsr = j.at("epsr");
                        lB = pc::lB( epsr );
                        depsdt = j.value("depsdt", -0.368*pc::T()/epsr);

                        ak.setRange(0, 1);
                        ak.setTolerance( j.value("tab_utol",1e-9) , j.value("tab_ftol",1e-2) );
                        bk.setRange(0, 1);
                        bk.setTolerance( j.value("tab_utol",1e-9) , j.value("tab_ftol",1e-2) );

                        if (type=="reactionfield") sfReactionField(j.at("eps_rf"));
                        if (type=="fanourgakis") sfFanourgakis();
                        if (type=="qpotential") sfQpotential(j.value("order",300));
			if (type=="q2potential") sfQ2potential(j.value("order",300));
                        if (type=="fennel") sfFennel(j.at("alpha"));
                        if (type=="plain") sfPlain(1);
                        if (type=="none") sfPlain(0);
                        if (type=="wolf") sfWolf(j.at("alpha"));

                        if ( tableA.empty() )
                            throw std::runtime_error("unknown coulomb type '" + type + "'" );
                    }

                    catch ( std::exception &e )
                    {
                        std::cerr << "DipoleDipoleGalore error: " << e.what();
                        throw;
                    }
                }
                
                DipoleDipoleGalore(double temperature, double epsr_in, double cutoff, string name, double parameter=0.0, double tab_utol=1e-9, double tab_ftol=1e-2) {
		  
                    try {
                        type = name;
                        name = "DipoleDipole-" + textio::toupper_first( type );
                        rc = cutoff;
                        rc2 = rc*rc;
                        rc1i = 1/rc;
                        epsr = epsr_in;
                        lB = pc::lB( epsr );
			pc::setT(temperature);

                        ak.setRange(0, 1);
                        ak.setTolerance( tab_utol , tab_ftol );
                        bk.setRange(0, 1);
                        bk.setTolerance( tab_utol , tab_ftol );

                        if (type=="reactionfield") sfReactionField(parameter);
                        if (type=="fanourgakis") sfFanourgakis();
                        if (type=="qpotential") sfQpotential(int(parameter));
			if (type=="q2potential") sfQ2potential(int(parameter));
                        if (type=="fennel") sfFennel(parameter);
                        if (type=="plain") sfPlain(1);
                        if (type=="none") sfPlain(0);
                        if (type=="wolf") sfWolf(parameter);

                        if ( tableA.empty() )
                            throw std::runtime_error("unknown coulomb type '" + type + "'" );
                    }

                    catch ( std::exception &e )
                    {
                        std::cerr << "DipoleDipoleGalore error: " << e.what();
                        throw;
                    }
                }

                template<class Tparticle>
                    double operator()(const Tparticle &a, const Tparticle &b, const Point &r) const {
                        double r1 = r.norm();
                        if (r1 < rc) {
                            double af = ak.eval(tableA,r1*rc1i);
                            double bf = bk.eval(tableB,r1*rc1i);
                            return lB*mu2mu(a.mu(), b.mu(), a.muscalar()*b.muscalar(), r,af,bf);
                        }
                        return 0.0;
                    }
                    
                /**
		 * @brief Self-energy of the potential
		 */
		template<class Tpvec>
		    double internal(const Tpvec &p, const Group &g) const { 
		    double Emu = 0;
		    for (auto i : g)
			Emu += p[i].mu().dot( p[i].mu() ) * p[i].muscalar() * p[i].muscalar();
		    return -selfenergy_prefactor*Emu*lB/rc/rc2;
		}

                double dielectric_constant(double M2V) {
                    return calcDielectric( M2V );
                } 

                string info(char w) {
                    using namespace textio;
                    std::ostringstream o;
                    o << pad(SUB, w, "Temperature") << pc::T() << " K" << endl
                        << pad(SUB, w, "Dielectric constant") << epsr << endl
                        << pad(SUB, w + 6, "T" + partial + epsilon + "/" + epsilon + partial + "T") << depsdt << endl
                        << pad(SUB, w, "Bjerrum length") << lB << _angstrom << endl
                        << pad(SUB,w,"Cutoff") << rc << _angstrom << endl;
                    if (type=="reactionfield") {
                        if(epsrf > 1e10)
                            o << pad(SUB,w+1, epsilon_m+"_RF") << infinity << endl;
                        else
                            o << pad(SUB,w+1, epsilon_m+"_RF") << epsrf << endl;
                    }
                    if (type=="qpotential" || type=="q2potential")
                        o << pad(SUB,w, "order") << order << endl;
                    if (type=="fennel" || type=="wolf")
                        o << pad(SUB,w, "alpha") << alpha << endl;
                    return o.str();
                }
        };

        /**
         * @brief Ion-dipole interaction energy
         *
         * \f$ u = q\boldsymbol{{\rm T}}_1(\boldsymbol{r})\cdot \boldsymbol{\mu} \f$ 
         * 
         * where
         * 
         * \f$ \boldsymbol{{\rm T}}_1(\boldsymbol{r}) = \nabla\left(\frac{1}{|\boldsymbol{r}|}\right). \f$
         * 
         * Here \f$ \boldsymbol{\mu} \f$ is the dipole moment, \f$ q \f$ is the charge and \f$ \boldsymbol{r} \f$ is the 
         * distance-vector between the two. The expression can also be written as
         * 
         * \f$ u = \frac {q\mu \cos(\theta)}{|\boldsymbol{r}|^2}  \f$
         * 
         * where \f$ \theta \f$ is the angle between \f$ \boldsymbol{r} \f$ and \f$ \boldsymbol{\mu} \f$.
         * 
         */
        class IonDipole : public PairPotentialBase {
            private:
                string _brief() { return "Ion-dipole"; }
            protected:
                double _lB;
            public:
                IonDipole(Tmjson &j) {
                    name="Ion-dipole";
                    _lB = Coulomb(j).bjerrumLength();
                }

                /**
                 * @brief Ion-dipole
                 * @param a Dipole particle A
                 * @param b Dipole particle B
                 * @param r Direction \f$ r_A - r_B \f$  
                 */
                template<class Tparticle>
                    double operator()(const Tparticle &a, const Tparticle &b, const Point &r) const {
                        return _lB*q2mu(a.charge*b.muscalar(),b.mu(),b.charge*a.muscalar(),a.mu(),r);
                    }

                string info(char w) { return _brief(); }
        };

        /**
         * @brief Ion-quadupole interaction energy
         *
         * \f$ u = q\boldsymbol{{\rm T}}_2(\boldsymbol{r}):\Theta \f$ 
         * 
         * where
         * 
         * \f$ \boldsymbol{{\rm T}}_2(\boldsymbol{r}) = \nabla^2\left(\frac{1}{|\boldsymbol{r}|}\right). \f$
         * 
         * Here \f$ q \f$ is the charge, \f$ \Theta \f$ the quadrupole moment and \f$ \boldsymbol{r} \f$ the 
         * distance-vector between the two. The expression can also be written as
         * 
         * \f$ u = 3\frac{\boldsymbol{r}^T\Theta \boldsymbol{r}}{|\boldsymbol{r}|^5} - \frac{tr(\Theta)}{|\boldsymbol{r}|^3}  \f$
         * 
         */
        class IonQuad : public PairPotentialBase {
            private:
                string _brief() { return "Ion-quadrupole"; }
            protected:
                double _lB;
            public:
                IonQuad(Tmjson &in) {
                    name="Ion-Quad";
                    _lB = Coulomb(in).bjerrumLength();
                }
                template<class Tparticle>
                    double operator()(const Tparticle &a, const Tparticle &b, const Point &r) const {
                        return _lB*q2quad(a.charge, b.theta,b.charge, a.theta,r);
                    }

                template<class Tparticle>
                    Point field(const Tparticle &p, const Point &r) const {
                        return Point(0,0,0);
                    }

                string info(char w) { return _brief(); }
        };

        /**
         * @brief Wolf summation for electrostatic interactions
         *
         * This will take care of long-ranged electrostatic interactions using the Wolf summation scheme. Currently it is possible to include ions, dipoles and quadrupoles. 
         * 
         *  Keyword          |  Description
         * :--------------   | :---------------
         * `kappa`           |  Damping-parameter                                                          
         * `cutoff`          |  Cut-off                                                                        
         * `forceshifted`    |  Sets if the potential will be force-shifted 
         * `tab_utol`        |  Tolerance for splines in energy-calculations.                         (Default: \f$ 10^{-9}\f$)
         * `tab_ftol`        |  Tolerance for splines in force-calculations.                          (Default: \f$ 10^{-2}\f$)
         * 
         * The formalism used the following basic interaction-tensor.
         * 
         * @f[
         * \boldsymbol{{\rm T}}_0(\boldsymbol{r}) = \frac{{\rm erfc}(\kappa |\boldsymbol{r}|)}{|\boldsymbol{r}|}
         * @f]
         * 
         * Higher order interaction-tensors are retrieved by taking the gradient of the previous one. I.e.
         * 
         * @f[
         * \boldsymbol{{\rm T}}_1(\boldsymbol{r}) = \nabla \boldsymbol{{\rm T}}_0(\boldsymbol{r}).
         * @f]
         * 
         * The effective potential is shifted and possibly also force-shifted according to the equation below.
         * 
         * @f[
         * \boldsymbol{{\rm T}}_{Mod}(\boldsymbol{r}) = \boldsymbol{{\rm T}}(\boldsymbol{r}) - \left.\boldsymbol{{\rm T}}(\boldsymbol{r})\right|_{|\boldsymbol{r}|\to R_c^-} - (|\boldsymbol{r}| - R_c)\frac{\partial \boldsymbol{{\rm T}}(\boldsymbol{r})}{\partial |\boldsymbol{r}|}
         * @f]
         * 
         * Original article for ionic interactions [doi](http://dx.doi.org/10.1063/1.478738)
         * Force-shifted exansion for ionic interactions [doi](http://dx.doi.org/10.1063/1.2206581)
         * Force-shifted exansion for dipolar interactions [doi](http://dx.doi.org/10.1063/1.4923001)
         * 
         * @todo Implement field-functions.
         */
        template<bool useIonIon=false, bool useIonDipole=false, bool useDipoleDipole=false, bool useIonQuadrupole=false>
            class MultipoleWolf : public PairPotentialBase {
                private:
                    double kappa, cutoff, cutoff2, der;
                    bool forceshifted;
                    Tabulate::Andrea<double> T0, T1, T2a, T2b; // splitting function
                    Tabulate::TabulatorBase<double>::data tableT0, tableT1, tableT2a, tableT2b; // data for splitting function

                    string _brief() {
                        std::ostringstream o;
                        o << "Multipole Wolf, lB=" << _lB << textio::_angstrom;
                        return o.str();          
                    }
                protected:
                    double _lB;

                public:
                    MultipoleWolf(Tmjson &in) {
                        name="Multipole Wolf";
                        _lB = Coulomb(in).bjerrumLength();

                        kappa = in.at("kappa");
                        cutoff = in.at("cutoff");
                        cutoff2 = cutoff*cutoff;
                        forceshifted = in.at("forceshifted");
                        der = 0.0;
                        if(forceshifted)
                            der = 1.0;
                        double kR = kappa*cutoff;

                        T0.setRange(0, 1);
                        T0.setTolerance(in.value("tab_utol",1e-9),in.value("tab_ftol",1e-2) );
                        tableT0 = T0.generate( [&](double q) { return ( erfc(kR*q) - q*erfc(kR) + der*q*(q - 1.0)*(erfc(kR)+2.0*kR/sqrt(pc::pi)*exp(-kR*kR)) ); } );

                        T1.setRange(0, 1);
                        T1.setTolerance(in.value("tab_utol",1e-9),in.value("tab_ftol",1e-2) );
                        tableT1 = T1.generate( [&](double q) { return ((2.0*kR*q/sqrt(pc::pi)*exp(-kR*kR*q*q) + erfc(kR*q)) - (2.0*kR/sqrt(pc::pi)*exp(-kR*kR) + erfc(kR))*q*q + der*q*q*(q - 1.0)*2.0*( (kR/sqrt(pc::pi)*exp(-kR*kR)*(1.0/sqrt(pc::pi) + 1.0 + 2.0*kR*kR/sqrt(pc::pi)) + erfc(kR) ))); } );

                        T2a.setRange(0, 1);
                        T2a.setTolerance(in.value("tab_utol",1e-9),in.value("tab_ftol",1e-2) );
                        tableT2a = T2a.generate( [&](double q) { return ( ( 2.0*kR*q*exp(-kR*kR*q*q)*(2.0*kR*kR*q*q/3.0 + 1.0)/sqrt(pc::pi) + erfc(kR*q)) - (erfc(kR) + 2.0*kR*exp(-kR*kR)*(2.0*kR*kR/3.0 + 1.0)/sqrt(pc::pi))*q*q*q + der*q*q*q*(q-1.0)*( 3.0*erfc(kR) + 2.0*kR*exp(-kR*kR)*(3.0 + 2.0*kR*kR + 4.0/3.0*kR*kR*kR*kR)/sqrt(pc::pi))); } );
			
                        T2b.setRange(0, 1);
                        T2b.setTolerance(in.value("tab_utol",1e-9),in.value("tab_ftol",1e-2) );
			tableT2b = T2b.generate( [&](double q) { return ( 4.0/3.0*kR*kR*kR*q*q*q*exp(-kR*kR*q*q)/sqrt(pc::pi) - 4.0/3.0*kR*kR*kR*exp(-kR*kR)/sqrt(pc::pi)*q*q*q + der*q*q*q*(q - 1.0)*8.0/3.0*exp(-kR*kR)*pow(kR,5)/sqrt(pc::pi) ); } );
		    }   

                    template<class Tparticle>
                        double operator()(const Tparticle &a, const Tparticle &b, const Point &r) {
                            double r2 = r.squaredNorm();
                            if(r2 > cutoff2)
                                return 0.0;

                            double U_total = 0;
                            double r1 = sqrt(r2);

                            if(useIonIon == true) U_total += a.charge * b.charge / r1 * T0.eval( tableT0, r1/cutoff );
                            if(useIonDipole == true) {
                                double T1_temp = T1.eval( tableT1, r1/cutoff );
                                U_total += a.charge*b.muscalar()*b.mu().dot(r)/r1*T1_temp/r2;
                                U_total += b.charge*a.muscalar()*a.mu().dot(-r)/r1*T1_temp/r2;
                            }
                            if(useDipoleDipole == true) U_total += mu2mu(a.mu(), b.mu(), a.muscalar()*b.muscalar(), r,T2a.eval(tableT2a,r1/cutoff),T2b.eval(tableT2b,r1/cutoff));
                            if(useIonQuadrupole == true) {
                                double traceA = a.theta.trace();
                                double traceB = b.theta.trace();
                                double crossA = r.transpose()*a.theta*r;
                                double crossB = r.transpose()*b.theta*r;
                                U_total += b.charge*((3*crossA/r2 - traceA)*T2a.eval(tableT2a,r1/cutoff) - traceA*T2b.eval(tableT2b,r1/cutoff))/r1/r2;
                                U_total += a.charge*((3*crossB/r2 - traceB)*T2a.eval(tableT2a,r1/cutoff) - traceB*T2b.eval(tableT2b,r1/cutoff))/r1/r2;
                            }
                            return _lB*U_total;
                        }

                    template<class Tparticle>
                        Point field(const Tparticle &p, const Point &r) {
                            return Point(0,0,0);
                        }

                    template<class Tparticle>
                        double fieldEnergy(const Tparticle &p, const Point &E) {
                            return 0;
                        }

                    string info(char w) {
                        using namespace textio;
                        std::ostringstream o;
                        o << pad(SUB,w,"Temperature") << pc::T() << " K" << endl
                            << pad(SUB,w,"Bjerrum length") << _lB << " "+angstrom << endl
                            << pad(SUB,w,"Cutoff") << cutoff << " "+angstrom << endl
                            << pad(SUB,w,"Kappa") << kappa << " "+angstrom+"^-1" << endl
                            << pad(SUB,w,"Force-shifted") << std::boolalpha << forceshifted << endl;
                        return o.str();
                    }
            };

        /**
         * @brief Ion-dipole interaction with long-ranged compensation using radial-derivates 1-3 to be zero at the cut-off for the original ionic interaction-tensor.
         * The potential also mimics the interaction-tensor \f$ T_0(r,\alpha) = erfc(\alpha r)/r \f$ with \f$ \alpha = \sqrt{\pi} \f$.
         * See DOI: 10.1021/jp510612w for more info.
         * 
         *  Keyword          |  Description
         * :--------------   | :---------------
         * `cutoff`          |  Cut-off for interactions.                               (Default: Infinity)
         * `eps_r`           |  Dielectric constant of the medium.                      (Default: \f$ \varepsilon_r = 1 \f$)
         */
        class IonDipoleSP3 : public IonDipole {
            private:
                string _brief() { return "IonDipole SP3"; }
                double rc1, rc1i, rc2;
            public:
                IonDipoleSP3(Tmjson &j) : IonDipole(j) { 
                    name += " SP3"; 
                    rc1 = j.at("cutoff");
                    rc1i = 1.0/rc1;
                    rc2 = rc1*rc1;
                }

                template<class Tparticle>
                    double operator()(const Tparticle &a, const Tparticle &b, const Point &r) const {
                        double r2 = r.squaredNorm();
                        if(r2 < rc2) {
                            double q = sqrt(r2)/rc1;
                            double q2 = q*q;
                            return IonDipole::operator()(a,b,r)*(1.0 - (21.0 - 35.0*q + 15.0*q2)*q2*q2*q);
                        }
                        return 0.0;
                    }

                string info(char w) {
                    using namespace textio;
                    std::ostringstream o;
                    o << IonDipole::info(w)
                        << pad(SUB,w,"Cutoff") << rc1 << " "+angstrom << endl;
                    return o.str();
                }
        };
    }
}
#endif

