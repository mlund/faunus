#ifndef FAUNUS_MULTIPOLE_H
#define FAUNUS_MULTIPOLE_H

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
            double ss = 1 - pow(2.71828,-std::min(expmax,pow((1/(asw*r1i)),nsw)));

            double uexp  = vec[0]*pow(2.71828,-std::min(expmax,vec[1]/r1i));
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
            double ss = 1 - pow(2.71828,-std::min(expmax,pow((1/(asw*r1i)),nsw)));

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
            double udis2  = vec[0]*pow(2.71828,-vec[1]/r1i);
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

            double uexp1  = vec[4]*pow(2.71828,-std::min(expmax,vec[5]/r1i));
            double uexp2  = 0;
            if(vec[6] != 0) uexp2 = vec[3]*pow(r1i,vec[6]);               // vec[6] = nab   <------------------------
            double udis1  =-vec[2]*r6i;
            double udis2  = vec[0]*pow(2.71828,-std::min(expmax,vec[1]/r1i));
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

            double uexp1  = vec[4]*pow(2.71828,-std::min(expmax,vec[5]/r1i));
            double uexp2  = 0;
            if(vec[6] != 0) uexp2 = vec[3]*pow(r1i,vec[6]);       // vec[6] = nab   <------------------------
            double udis1  =-vec[2]*r6i;
            double udd = 1 + ud1 + ud2 + ud3 + ud4 + ud5 + ud6;
            double udis2  = vec[0]*pow(2.71828,-std::min(expmax,1/bri));
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

            double uexp1  = vec[4]*pow(2.71828,-std::min(expmax,vec[5]/r1i));
            double uexp2  = 0;
            if(vec[6] != 0) uexp2 = vec[3]*pow(r1i,vec[6]);                     // vec[6] = nab   <------------------------
            double udis1  =-vec[2]*r6i;
            double udd = 1 + ud1 + ud2 + ud3 + ud4 + ud5 + ud6;
            double udis2  = vec[0]*pow(2.71828,-std::min(expmax,1/bri));
            double uchtexp = -vec[8]*pow(2.71828,-std::min(expmax,vec[7]/r1i));  // vec[7] = acht, vec[8] = kcht    <------------------------
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

            double uchtexp = -vec[8]*pow(2.71828,-std::min(expmax,vec[7]*(pow((r.norm()-vec[3]),2))));  // vec[7] = acht, vec[8] = kcht    <------------------------
            double uexp  = vec[4]*pow(2.71828,-std::min(expmax,vec[5]/r1i));
            double udis1  =-vec[2]*r6i;
            double udd = 1 + ud1 + ud2 + ud3 + ud4 + ud5 + ud6;
            double udis2  = vec[0]*pow(2.71828,-std::min(expmax,1/bri));
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
        double mu2mu(const Tvec &muA, const Tvec &muB, double muAxmuB, const Tvec &r) {
#ifdef FAU_APPROXMATH
            double r1i = invsqrtQuake( r.squaredNorm() );
            double r2i = r1i*r1i;
#else
            double r2i = 1/r.squaredNorm();
            double r1i = sqrt(r2i);
#endif
            double r3i = r1i*r2i;
            //double r5i = r3i*r2i;
            //Eigen::Matrix3d T = 3*r5i*r*r.transpose() - r3i*Matrix3d::Identity();
            //double W = -muA.transpose()*T*muB;
            double W = r3i*(3*muA.dot(r)*muB.dot(r)*r2i - muA.dot(muB));
            return -W*muAxmuB;
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

    /**
     * @brief Base class for Wolf based interactions
     *
     * The idea is that this class has no dependencies and is
     * to be used as a helper class for other classes.
     */
    class WolfBase {

        struct wdata {
            double r1i_d, r2i, der_dT0c, T1, T1c_r1i, der_dT1c_r1i, T21, T22, T2c2_r2i, der_dT2c1, der_dT2c2_r2i;
        };
        private:
        double rc1, rc1i, rc1i_d, rc2i, kappa, kappa2, constant;
        double dT0c, T1c_rc1, dT1c_rc1, T2c1, T2c2_rc2, dT2c1, dT2c2_rc2;
        wdata data;

        public:
        void calcWolfData(const Point &r) {
            data.r2i = 1/r.squaredNorm();
            double r1i = sqrt(data.r2i);
            data.r1i_d = erfc_x(kappa/r1i)*r1i;
            double der = (1/r1i) - rc1;
            der = 0;
            double expK = constant*exp(-kappa2/data.r2i);
            data.der_dT0c = der*dT0c;
            data.T1 = (expK + data.r1i_d)*data.r2i;
            data.T1c_r1i = T1c_rc1*r1i;
            data.der_dT1c_r1i = der*dT1c_rc1*r1i;
            data.T21 = -(data.r1i_d + expK)*data.r2i;
            data.T22 = (3.*data.r1i_d*data.r2i + (3.*data.r2i + 2.*kappa2)*expK)*data.r2i;
            data.T2c2_r2i = T2c2_rc2*data.r2i;
            data.der_dT2c1 = der*dT2c1;
            data.der_dT2c2_r2i = der*dT2c2_rc2*data.r2i;
            data.T2c2_r2i = 0.0; // No Energy
            data.der_dT2c2_r2i = 0.0; // No Force
        }

        /**
         * @brief Constructor
         * @param alpha Dampening factor (inverse angstrom)
         * @param rcut Cutoff distance (angstrom)
         */
        WolfBase(double alpha, double rcut) {
            kappa = alpha;
            kappa2 = kappa*kappa;
            constant = 2*kappa/sqrt(pc::pi);
            rc1 = rcut;
            double rc2 = rc1*rc1;
            rc2i = 1/rc2;
            rc1i = 1/rc1;
            double expKc = constant*exp(-kappa2/rc2i);
            rc1i_d = erfc_x(kappa*rc1)*rc1i;

            T1c_rc1 = (expKc + rc1i_d)*rc2i;
            T2c1 = -(expKc + rc1i_d)*rc2i;
            T2c2_rc2 = (3.*rc1i_d*rc2i*rc2i + (3.*rc2i + 2.*kappa2)*expKc*rc2i);

            dT0c = -(expKc + rc1i_d)*rc1i;
            dT1c_rc1 = (-2*T1c_rc1/rc1) - 2*kappa2*expKc*rc1i;
            dT2c1 = -(3*T2c1/rc1) + (2*kappa2*exp(-rc2*kappa2)*rc1i*constant);
            dT2c2_rc2 = (-3*T2c2_rc2/rc1) - (4*kappa2*kappa2*exp(-rc2*kappa2)*rc1i*constant);

            T1c_rc1 = T1c_rc1*rc1;
            dT1c_rc1 = dT1c_rc1*rc1;
            T2c2_rc2 = T2c2_rc2*rc2;
            dT2c2_rc2 = dT2c2_rc2*rc2;
        }

        /**
         * @brief Returns ion-ion interaction.
         * @param qA Charge of ion A
         * @param qB Charge of ion B
         * @param r Direction \f$ r_A - r_B \f$
         */
        template<bool useWdata=false, class Tvec>
            double q2q(double qA, double qB, const Tvec &r) const {
                /* Code to use if calcWolfData is not called */
                if (useWdata==false) {
                    double r2i = 1/r.squaredNorm();
                    if (r2i < rc2i)
                        return 0;
                    double r1i = sqrt(r2i);
                    double der = (1/r1i) - rc1;
                    double r1i_d = erfc_x(kappa/r1i)*r1i;
                    return (qA*qB*(r1i_d - rc1i_d - der*dT0c));
                }
                return (qA*qB*(data.r1i_d - rc1i_d - data.der_dT0c));
            }

        /**
         * @brief Returns ion-dipole interaction.
         * @param QBxMuA Product of ion B:s charge and dipole A:s scalar
         * @param muA Unit dipole moment vector of particel A
         * @param QAxMuB Product of ion A:s charge and dipole B:s scalar
         * @param muB Unit dipole moment vector of particel B
         * @param r Direction \f$ r_A - r_B \f$
         */
        template<bool useWdata=false, class Tvec>
            double q2mu(double QBxMuA, const Tvec &muA, double QAxMuB, const Tvec &muB, const Tvec &r) const {
                /* Code to use if calcWolfData is not called */
                if (useWdata==false) {
                    double r2i = 1/r.squaredNorm();
                    if (r2i < rc2i)
                        return 0;
                    double r1i = sqrt(r2i);
                    double T1 = (constant*exp(-kappa2/r2i) + erfc_x(kappa/r1i)*r1i)*r2i;
                    double der = (1/r1i) - rc1;
                    double W1 = QBxMuA*muA.dot(r)*(T1 - T1c_rc1*r1i - der*dT1c_rc1*r1i);
                    double W2 = QAxMuB*muB.dot(-r)*(T1 - T1c_rc1*r1i - der*dT1c_rc1*r1i);
                    return (W1 + W2);
                }
                double W1 = QBxMuA*muA.dot(r)*(data.T1 - data.T1c_r1i - data.der_dT1c_r1i);
                double W2 = QAxMuB*muB.dot(-r)*(data.T1 - data.T1c_r1i - data.der_dT1c_r1i);
                return (W1 + W2);
            }

        /**
         * @brief Dipole-dipole energy
         * @param muA Dipole moment (A) unit vector
         * @param muB Dipole moment (B) unit vector
         * @param muAxmuB Product of dipole moment scalars, |A|*|B|
         * @param r Direction \f$ r_A - r_B \f$
         */
        template<bool useWdata=false, class Tvec>
            double mu2mu(const Tvec &muA, const Tvec &muB, double muAxmuB, const Tvec &r) const {
                /* Code to use if calcWolfData is not called */
                if (useWdata==false) {
                    double r2i = 1.0/r.squaredNorm();
                    if (r2i < rc2i)
                        return 0;
                    double r1i = sqrt(r2i);
                    double r1i_d = erfc_x(kappa/r1i)*r1i;
                    double expK = constant*exp(-kappa2/r2i);
                    double T2_1 = -(r1i_d + expK)*r2i;
                    double T2_2 = (3.*r1i_d*r2i + (3.*r2i + 2.*kappa2)*expK)*r2i;
                    double der = (1/r1i) - rc1;
                    double t3 = muA.dot(muB)*(T2_1 - T2c1 - der*dT2c1);
                    double t5 = muA.dot(r)*muB.dot(r)*(T2_2 - T2c2_rc2*r2i - der*dT2c2_rc2*r2i);
                    return -(t5 + t3)*muAxmuB;
                }
                double t3 = muA.dot(muB)*(data.T21 - T2c1 - data.der_dT2c1);
                double t5 = muA.dot(r)*muB.dot(r)*(data.T22 - data.T2c2_r2i - data.der_dT2c2_r2i);
                return -(t5 + t3)*muAxmuB;
            }

        /**
         * @brief Returns ion-quadrupole energy
         * @param qA Charge of particle A
         * @param quadB Quadrupole moment of particle B
         * @param qB Charge of particle B
         * @param quadA Quadrupole moment of particle A
         * @param r Direction @f$ r_A - r_B @f$
         */
        template<bool useWdata=false, class Tvec, class Tmat>
            double q2quad(double qA, const Tmat &quadB,double qB, const Tmat &quadA, const Tvec &r) const {
                /* Code to use if calcWolfData is not called */
                if (useWdata==false) {
                    double r2i = 1/r.squaredNorm();
                    if (r2i < rc2i)
                        return 0;
                    double r1i = sqrt(r2i);
                    double r1i_d = erfc_x(kappa/r1i)*r1i;
                    double expK = constant*exp(-kappa2/r2i);
                    double T2_1 = -(r1i_d + expK)*r2i;
                    double T2_2 = (3.*r1i_d*r2i + (3.*r2i + 2.*kappa2)*expK)*r2i;
                    double der = (1/r1i) - rc1;
                    double WAB = r.transpose()*quadB*r;
                    WAB = WAB*(T2_2 - T2c2_rc2*r2i - der*dT2c2_rc2*r2i) + quadB.trace()*(T2_1 - T2c1 - der*dT2c1);
                    double WBA = r.transpose()*quadA*r;
                    WBA = WBA*(T2_2 - T2c2_rc2*r2i - der*dT2c2_rc2*r2i) + quadA.trace()*(T2_1 - T2c1 - der*dT2c1);
                    return (qA*WAB + qB*WBA);
                }
                double WAB = r.transpose()*quadB*r;
                WAB = WAB*(data.T22 - data.T2c2_r2i - data.der_dT2c2_r2i) + quadB.trace()*(data.T21 - T2c1 - data.der_dT2c1);
                double WBA = r.transpose()*quadA*r;
                WBA = WBA*(data.T22 - data.T2c2_r2i - data.der_dT2c2_r2i) + quadA.trace()*(data.T21 - T2c1 - data.der_dT2c1);
                return (qA*WAB + qB*WBA);
            }

        /** 
         * @brief Field at `r` due to charge `p` 
         * @param p Particles from which field arises
         * @param r Direction @f$ r_A - r_B @f$
         */
        template<bool useWdata=false, class Tparticle>
            Point fieldCharge(const Tparticle &p, const Point &r) const {
                /* Code to use if calcWolfData is not called */
                if (useWdata==false) {
                    double r2i = 1/r.squaredNorm();
                    if (r2i < rc2i)
                        return Point(0,0,0);
                    double r1i = sqrt(r2i);
                    double T1 = (constant*exp(-kappa2/r2i) + erfc_x(kappa/r1i)*r1i)*r2i;
                    double der = (1/r1i) - rc1;
                    return (T1 - T1c_rc1*r1i - der*dT1c_rc1*r1i)*r*p.charge;
                }
                return (data.T1 - data.T1c_r1i - data.der_dT1c_r1i)*r*p.charge;
            }

        /** 
         * @brief Field at `r` due to dipole `p` 
         * @param p Particles from which field arises
         * @param r Direction @f$ r_A - r_B @f$
         */
        template<bool useWdata=false, class Tparticle>
            Point fieldDipole(const Tparticle &p, const Point &r) const {
                /* Code to use if calcWolfData is not called */
                if (useWdata==false) {
                    double r2i = 1.0/r.squaredNorm();
                    if (r2i < rc2i)
                        return Point(0,0,0);
                    double r1i = sqrt(r2i);
                    double r1i_d = erfc_x(kappa/r1i)*r1i;
                    double expK = constant*exp(-kappa2/r2i);
                    double T2_1 = -(r1i_d + expK)*r2i;
                    double T2_2 = (3.*r1i_d*r2i + (3.*r2i + 2.*kappa2)*expK)*r2i;
                    double der = (1/r1i) - rc1;
                    Point t3 = p.mu*(T2_1 - T2c1 - der*dT2c1);
                    Point t5 = r*p.mu.dot(r)*(T2_2 - T2c2_rc2*r2i - der*dT2c2_rc2*r2i);
                    return (t5 + t3)*p.muscalar;
                }
                Point t3 = p.mu*(data.T21 - T2c1 - data.der_dT2c1);
                Point t5 = r*p.mu.dot(r)*(data.T22 - data.T2c2_r2i - data.der_dT2c2_r2i);
                return (t5 + t3)*p.muscalar;
            }

        double getRc2i() const { return rc2i; }
        double getR2i() const { return data.r2i; }
        double getKappa() { return kappa; }
        double getCutoff() { return rc1; }
    };


    /**
     * @brief Base class for Gaussian-damped interactions. Implemented according to DOI: 10.1002/jcc.20574
     *
     * The idea is that this class has no dependencies and is
     * to be used as a helper class for other classes.
     */
    class GaussianDampingBase {
        private:
            double constant;
            Eigen::VectorXd betaC, betaD, betaQ, betaC3, betaD2, betaD3;
            Eigen::MatrixXd betaCC, betaCD, betaCQ, betaDD;
            Eigen::MatrixXd betaCC2, betaCD2, betaCQ2, betaDD2;
            Eigen::MatrixXd betaCC3, betaCD3, betaCQ3, betaDD3;

        public:
            /**
             * @brief Constructor
             */
            GaussianDampingBase() {
                constant = 2/sqrt(pc::pi);
                int N = atom.size() - 1;

                double alpha;
                double pre_factor = pow(3*sqrt(8*pc::pi)/4,1.0/3.0);
                for(int i = 0; i < N; i++) {
                    alpha = (atom[i+1].alpha(0,0) + atom[i+1].alpha(1,1) + atom[i+1].alpha(2,2))/3.0;
                    if(atom[i+1].betaC == pc::infty) {
                        atom[i+1].betaC = 0.75*pre_factor*pow(alpha,-1.0/3.0);
                    }
                    if(atom[i+1].betaD == pc::infty) {
                        atom[i+1].betaD = 0.75*pre_factor*pow(alpha,-1.0/3.0);
                    }
                    if(atom[i+1].betaQ == pc::infty) {
                        atom[i+1].betaC = 0.75*pre_factor*pow(alpha,-1.0/3.0);
                    }
                }

                betaC.resize(N);
                betaD.resize(N);
                betaQ.resize(N);
                betaC3.resize(N);
                betaD2.resize(N);
                betaD3.resize(N);
                betaCC.resize(N,N);
                betaCD.resize(N,N);
                betaCQ.resize(N,N);
                betaDD.resize(N,N);
                betaCC2.resize(N,N);
                betaCD2.resize(N,N);
                betaCQ2.resize(N,N);
                betaDD2.resize(N,N);
                betaCC3.resize(N,N);
                betaCD3.resize(N,N);
                betaCQ3.resize(N,N);
                betaDD3.resize(N,N);

                for(int i = 0; i < N; i++) {
                    betaC(i) = atom[i+1].betaC;
                    betaC3(i) = betaC(i)*betaC(i)*betaC(i);
                    betaD(i) = atom[i+1].betaD;
                    betaD2(i) = betaD(i)*betaD(i);
                    betaD3(i) = betaD2(i)*betaD(i);
                    betaQ(i) = atom[i+1].betaQ;
                }
                for(int i = 0; i < N; i++) {
                    for(int j = i; j < N; j++) {
                        betaCC(i,j) = betaC(i)*betaC(j)/sqrt(betaC(i)*betaC(i) + betaC(j)*betaC(j));
                        betaCD(i,j) = betaC(i)*betaD(j)/sqrt(betaC(i)*betaC(i) + betaD(j)*betaD(j));
                        betaCQ(i,j) = betaC(i)*betaQ(j)/sqrt(betaC(i)*betaC(i) + betaQ(j)*betaQ(j));
                        betaDD(i,j) = betaD(i)*betaD(j)/sqrt(betaD(i)*betaD(i) + betaD(j)*betaD(j));
                        betaCC(j,i) = betaCC(i,j);
                        betaCD(j,i) = betaD(i)*betaC(j)/sqrt(betaD(i)*betaD(i) + betaC(j)*betaC(j));
                        betaCQ(j,i) = betaQ(i)*betaC(j)/sqrt(betaQ(i)*betaQ(i) + betaC(j)*betaC(j));
                        betaDD(j,i) = betaDD(i,j);
                        betaCC2(i,j) = betaCC(i,j)*betaCC(i,j);
                        betaCD2(i,j) = betaCD(i,j)*betaCD(i,j);
                        betaCQ2(i,j) = betaCQ(i,j)*betaCQ(i,j);
                        betaDD2(i,j) = betaDD(i,j)*betaDD(i,j);
                        betaCC2(j,i) = betaCC2(i,j);
                        betaCD2(j,i) = betaCD(j,i)*betaCD(j,i);
                        betaCQ2(j,i) = betaCQ(j,i)*betaCQ(j,i);
                        betaDD2(j,i) = betaDD2(i,j);
                        betaCC3(i,j) = betaCC2(i,j)*betaCC(i,j);
                        betaCD3(i,j) = betaCD2(i,j)*betaCD(i,j);
                        betaCQ3(i,j) = betaCQ2(i,j)*betaCQ(i,j);
                        betaDD3(i,j) = betaDD2(i,j)*betaDD(i,j);
                        betaCC3(j,i) = betaCC3(i,j);
                        betaCD3(j,i) = betaCD2(j,i)*betaCD(j,i);
                        betaCQ3(j,i) = betaCQ2(j,i)*betaCQ(j,i);
                        betaDD3(j,i) = betaDD3(i,j);
                    }
                }
            }

            /**
             * @brief Returns ion-ion energy
             * @param qA Charge of ion A
             * @param qB Charge of ion B
             * @param ida Id of particle A
             * @param idb Id of particle B
             * @param r Direction \f$ r_A - r_B \f$
             */
            template<class Tvec>
                double q2q(double qA, double qB, int ida, int idb, const Tvec &r) const {
                    double r1 = r.norm();
                    return (qA*qB*erf_x(betaCC(ida-1,idb-1)*r1)/r1);
                }

            /**
             * @brief Returns ion-dipole energy
             * @param QBxMuA Product of charge of particle B and dipole scalar of particle A
             * @param muA Unit dipole moment vector of particle A
             * @param QAxMuB Product of charge of particle A and dipole scalar of particle B
             * @param muB Unit dipole moment vector of particle B
             * @param ida Id of particle A
             * @param idb Id of particle B
             * @param r Direction \f$ r_A - r_B \f$
             */
            template<class Tvec>
                double q2mu(double QBxMuA, const Tvec &muA, double QAxMuB, const Tvec &muB, int ida, int idb, const Tvec &r) const {
                    double r2 = r.squaredNorm();
                    double r1 = sqrt(r2);
                    double B1_BA = ( (erf_x( betaCD(ida-1,idb-1) * r1 )/r1) - betaCD(ida-1,idb-1) * constant * exp(-betaCD2(ida-1,idb-1)*r2) ) / r2;
                    double B1_AB = ( (erf_x( betaCD(idb-1,ida-1) * r1 )/r1) - betaCD(idb-1,ida-1) * constant * exp(-betaCD2(idb-1,ida-1)*r2) ) / r2;
                    double W_BA = ( QBxMuA * muA.dot(r) * B1_BA);
                    double W_AB = ( QAxMuB * muB.dot(-r) * B1_AB);
                    return ( W_BA + W_AB );
                }

            /**
             * @brief Dipole-dipole energy
             * @param muA Unit dipole moment vector of particle A
             * @param muB Unit dipole moment vector of particle B
             * @param muAxmuB Product of dipole moment scalars
             * @param ida Id of particle A
             * @param idb Id of particle B
             * @param r Direction \f$ r_A - r_B \f$
             */
            template<class Tvec>
                double mu2mu(const Tvec &muA, const Tvec &muB, double muAxmuB, int ida, int idb, const Tvec &r) const {
                    double x = betaDD(ida-1,idb-1)*r.norm();
                    double x2 = x*x;
                    double erfX = erf_x(x)/x;
                    double expX = constant * exp(-x2);
                    double B1 = ( erfX - expX ) / x2;
                    double B2 = ( 3*erfX - (3 + 2*x2) * expX ) / ( x2 * x2 );
                    double W = ( muA.dot(muB) * B1 - betaDD2(ida-1,idb-1) * ( muA.dot(r) ) * ( muB.dot(r) ) * B2 ) * betaDD3(ida-1,idb-1);
                    return ( muAxmuB * W );
                }

            /**
             * @brief Returns ion-quadrupole energy
             * @param qA Charge of particle A
             * @param quadB Quadrupole moment of particle B
             * @param qB Charge of particle B
             * @param quadA Quadrupole moment of particle A
             * @param ida Id of particle A
             * @param idb Id of particle B
             * @param r Direction @f$ r_A - r_B @f$
             */
            template<class Tvec, class Tmat>
                double q2quad(double qA, const Tmat &quadB,double qB, const Tmat &quadA, int ida, int idb, const Tvec &r) const {
                    double r1 = r.norm();
                    double x_AB = betaCQ(ida-1,idb-1)*r1;
                    double x2_AB = x_AB*x_AB;
                    double erfX_AB = erf_x(x_AB)/x_AB;
                    double expX_AB = constant * exp(-x2_AB);
                    double B1_AB = ( erfX_AB - expX_AB ) / x2_AB;
                    double B2_AB = betaCQ2(ida-1,idb-1) * ( 3*erfX_AB - (3 + 2*x2_AB) * expX_AB ) / ( x2_AB * x2_AB );
                    double W_AB = r.transpose()*quadB*r;
                    W_AB = W_AB * B2_AB - quadB.trace() * B1_AB;
                    double x_BA = betaCQ(idb-1,ida-1)*r1;
                    double x2_BA = x_BA*x_BA;
                    double erfX_BA = erf_x(x_BA)/x_BA;
                    double expX_BA = constant * exp(-x2_BA);
                    double B1_BA = ( erfX_BA - expX_BA ) / x2_BA;
                    double B2_BA = betaCQ2(idb-1,ida-1) * ( 3*erfX_BA - (3 + 2*x2_BA) * expX_BA ) / ( x2_BA * x2_BA );
                    double W_BA = r.transpose()*quadA*r;
                    W_BA = W_BA * B2_BA - quadA.trace() * B1_BA;
                    return ( qA * W_AB * betaCQ3(ida-1,idb-1) + qB * W_BA * betaCQ3(idb-1,ida-1) );
                }

            /** 
             * @brief Field at `r` due to ion `p`.
             * @param p Particles from which field arises
             * @param r Direction @f$ r_A - r_B @f$
             * @param ida Id of particle exposed to the field from p. If ida is not set the exposed particle is assumed to be a point particle (optional)
             */
            template<class Tparticle>
                Point fieldCharge(const Tparticle &p, const Point &r, int ida=-1) const {
                    Point E(0,0,0);
                    if(ida != -1) {
                        double x = betaCC(ida-1,p.id-1)*r.norm();
                        double x2 = x*x;
                        E = (p.charge * betaCC3(ida-1,p.id-1) * ( ( erf_x(x) / x ) - constant * exp(-x2) ) / x2)*r;
                    } else {
                        double x = betaC(p.id-1)*r.norm();
                        double x2 = x*x;
                        E = (p.charge * betaC3(p.id-1)* ( ( erf_x(x) / x ) - constant * exp(-x2) ) / x2) *r;
                    }
                    return E;
                }

            /** 
             * @brief Field at `r` due to dipole `p`.
             * @param p Particles from which field arises
             * @param r Direction @f$ r_A - r_B @f$
             * @param ida Id of particle exposed to the field from p. If ida is not set the exposed particle is assumed to be a point particle (optional)
             */
            template<class Tparticle>
                Point fieldDipole(const Tparticle &p, const Point &r, int ida=-1) const {
                    Point E(0,0,0);
                    if(ida != -1) {
                        double x = betaDD(ida-1,p.id-1)*r.norm();
                        double x2 = x*x;
                        double erfX = erf_x(x)/x;
                        double expX = constant * exp(-x2);
                        double B1 = ( erfX - expX ) / x2;
                        double B2 = ( 3*erfX - (3 + 2*x2) * expX ) / ( x2 * x2 );
                        E = -p.muscalar*( B1 * p.mu - betaDD2(ida-1,p.id-1) * p.mu.dot(r) * B2 *r ) * betaDD3(ida-1,p.id-1);
                    } else {
                        double x = betaD(p.id-1)*r.norm();
                        double x2 = x*x;
                        double erfX = erf_x(x)/x;
                        double expX = constant * exp(-x2);
                        double B1 = ( erfX - expX ) / x2;
                        double B2 = ( 3*erfX - (3 + 2*x2) * expX ) / ( x2 * x2 );
                        E = -p.muscalar*( B1 * p.mu - betaD2(p.id-1) * p.mu.dot(r) * B2 *r ) * betaD3(p.id-1);
                    }
                    return E;
                }
    };

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
         * @brief Ion-dipole interaction, 
         *
         * More info...
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
         * @brief Dipole-dipole interaction
         *
         * More info...
         */
        class DipoleDipole : public PairPotentialBase {
            private:
                string _brief() {
                    std::ostringstream o;
                    o << "Dipole-dipole, lB=" << _lB << textio::_angstrom;
                    return o.str();          
                }
            protected:
                double _lB;
            public:
                DipoleDipole(double T_kelvin, double epsilon_r) {
                    pc::setT(T_kelvin);
                    _lB = pc::lB(epsilon_r);
                }

                DipoleDipole(Tmjson &j) {
                    name="Dipole-dipole";
                    _lB = Coulomb(j).bjerrumLength();
                }
                template<class Tparticle>
                    double operator()(const Tparticle &a, const Tparticle &b, const Point &r) const {
                        return _lB*mu2mu(a.mu(), b.mu(), a.muscalar()*b.muscalar(), r);
                    }

                /** @brief Dipole field at `r` due to dipole `p` 
                 */
                template<class Tparticle>
                    Point field(const Tparticle &p, const Point &r) const {
                        double r2i = 1.0/r.squaredNorm();
                        double r1i = sqrt(r2i);
                        return ((3.0*p.mu().dot(r)*r*r2i - p.mu())*r2i*r1i)*p.muscalar()*_lB;
                    }

                /**
                 * @brief Interaction of dipole `p` with field `E`, see 'Intermolecular and SUrface Forces' by J. Israelachvili, p. 97 eq. 5.15
                 * @todo Needs to be tested!
                 */
                template<class Tparticle>
                    double fieldEnergy(const Tparticle &p, const Point &E) {
                        return -p.muscalar()*p.mu().dot(E);
                    }

                string info(char w) {
                    using namespace textio;
                    std::ostringstream o;
                    o << pad(SUB,w,"Temperature") << pc::T() << " K" << endl
                        << pad(SUB,w,"Bjerrum length") << _lB << " "+angstrom << endl;
                    return o.str();
                }
        };

        /** @brief Ion-quadrupole interaction
	 * 
         * More info...
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
         * @brief Coulomb type potentials with spherical cutoff
         *
         * Beyond a spherical cutoff, \f$r_c\f$, the potential is cut and shifted to zero while
         * \f$ u(r) = \lambda_B z_i z_j \mathcal{S} / r \f$
         * is used otherwise with the following choice of splitting functions, \f$\mathcal{S}\f$,
         *
         *  Coulomb-type     |  Splitting function, \f$\mathcal{S}\f$ | Additional keywords  | Reference
         *  ---------------- |  ------------------------------------- | -------------------- | ---------
         *  `yukawa`         |  \f$\exp(-\kappa r) + ...\f$           | `debyelength`        | doi
         *  `reactionfield`  |  ...                                   | `epsrf`              | doi
         *  `qpotential`     |  ...                                   | ...                  | doi
         *  `yonezawa`       |  ...                                   | ...                  | http://dx.doi.org/10/j97
         *  `grekiss`        |  ...                                   | ...                  | doi
         *  `plain`          |  \f$ 1 + ... \f$                       | none                 | doi
         *
         *  The following keywords are *required*:
         *
         *  Keyword       |  Description
         *  ------------- |  -------------------------------------------
         *  `coulombtype` |  Type of splitting function as defined above
         *  `cutoff`      |  Spherical cutoff in angstroms
         *  `epsr`        |  Relative dielectric constant of the medium
         *
         *  @todo Under construction!
         */
        class CoulombGalore : public PairPotentialBase {
            private:
                Tabulate::Andrea<double> sf; // splitting function
                Tabulate::TabulatorBase<double>::data table; // data for splitting function

                std::function<double(double)> calcDielectric; // function for dielectric const. calc.

                double lB, rc, rc2, rci, epsr;

                void sfReactionField(const Tmjson &j) {
                    double epsrf = j.at("epsrf");
                    double A = ( epsrf - epsr ) / ( 2 * epsrf + epsr );
                    double B = -3 * epsrf / ( 2 * epsrf + epsr );
                    table = sf.generate( [&](double q) { return 1 + A*q*q*q + B*q ; } ); 
                    calcDielectric = [&](double M2V) {
                        return 0.5 * ( 2 * epsrf - 1
                                + sqrt( -72*M2V*M2V*epsrf + 4*epsrf*epsrf + 4*epsrf + 1)
                                ) / (3 * M2V - 1 ); // Needs to be checked!
                    };
                    // we could also fill in some info string or JSON output...
                }

                void sfYukawa(const Tmjson &j) {
                    double kappa = 1 / j.at("debyelength").get<double>();
                    table = sf.generate( [&](double q) { return std::exp(-q*rc*kappa) ; } ); 
                    // we could also fill in some info string or JSON output...
                }

            public:
                CoulombGalore(const Tmjson &j) {

                    string type = j.at("coulombtype");
                    rc = j.at("cutoff");
                    rc2 = rc*rc;
                    rci = 1/rc;
                    epsr = j.at("epsr");
                    lB = pc::lB( epsr );

                    sf.setRange(0, 1);
                    sf.setTolerance(
                            j.value("tab_utol",1e-9),j.value("tab_ftol",1e-2) );

                    if (type=="reactionfield")
                        sfReactionField(j);

                    if (type=="yukawa")
                        sfYukawa(j);
                }

                template<class Tparticle>
                    double operator()(const Tparticle &a, const Tparticle &b, double r2) {
                        if (r2 < rc2) {
                            double r = sqrt(r2);
                            return lB * a.charge * b.charge / r * sf.eval( table, r*rci );
                        }
                        return 0;
                    }

                template<class Tparticle>
                    double operator()(const Tparticle &a, const Tparticle &b, const Point &r) const {
                        return operator()(a,b,r.squaredNorm());
                    }
        };

        /**
         * @brief Ion-ion interaction with reaction field. The potential can be shifted such to be zero at the cut-off.
         * 
         *  Keyword          |  Description
         * :--------------   | :---------------
         * `cutoff`          |  Cut-off for interactions.
         * `epsr`            |  Dielectric constant of the medium.                      ( Default: 1 )
         * `eps_rf`          |  Dielectric constant of the surroundings.
         * `kappa`           |  Inverse Debye screening length.                         ( Default: 0 )
         * `tab_utol`        |  Tolerance in splitting-function error.                  (Default: \f$ 10^{-9}\f$)
         * `tab_ftol`        |  Tolerance of splitting-function derivative error.       (Default: \f$ 10^{-2}\f$)
         * 
         * @note If 'eps_rf' is set to ( epsr , < 0 , 0 ) then 'vacuum'/insulating/conducting boundary conditions will be used. 
         * @warning Untested! Generalized approach using ionic strength does not work!
         * 
	 * @DOI See DOI: 10.1080/00268978000100361 for specifics on the reaction field potential.
         * @DOI See DOI: 10.1080/00268978300102721 for specifics on the dielectric constant.
         * @DOI See DOI: 10.1063/1.469273 for the generalized approach using ionic strength in the surrounding.
         */
        class CoulombRF : public Coulomb {
            private:
                string _brief() { return "Coulomb (reaction field)"; }
                double rc1,rc2,eps_RF,epsr,krf,crf,kappa,kappa2, prefactorA, prefactorB;
                bool eps_inf, eps_ins, eps_vac, eps_user; // Surrounding is: Conducting, Insulating, 'Vacuum', Set by user

                Tabulate::Andrea<double> splitting_function;
                Tabulate::TabulatorBase<double>::data table;
            public:
                CoulombRF(Tmjson &j) : Coulomb(j) {
                    name+=" Reaction Field";
                    rc1 = j.at("cutoff");
                    rc2 = rc1*rc1;
                    epsr = j.value("epsr",1.0);
                    eps_RF = j.at("eps_rf");
                    kappa = j.value("kappa",0.0);
                    kappa2 = kappa*kappa;
                    eps_inf = false;
                    eps_ins = false;
                    eps_user = false;
                    eps_vac = false;
                    if(fabs(eps_RF) < 1e-6) {
                        eps_RF = pc::infty;
                        eps_inf = true; // Conducting boundary conditions
                        prefactorA = 0.5;
			prefactorB = -1.5;
                    } else if(fabs(eps_RF-epsr) < 1e-6) {
                        eps_ins = true; // Insulating boundary conditions
                        prefactorA = 0.0;
			prefactorB = -1.0;
                    } else if(fabs(eps_RF-1.0) < 1e-6) {
                        eps_vac = true; // Vacuum boundary conditions
                        prefactorA = (1.0 - epsr)/(2.0+epsr);
			prefactorB = -3.0/(2.0 + epsr);
                    } else {
                        eps_user = true; // Set by user
                        prefactorA = (eps_RF - epsr)/(2.0*eps_RF+epsr);
			prefactorB = -3.0*eps_RF/(2.0*eps_RF + epsr);
                    }

                    std::function<double(double)> PFF = [&](double q) { return (1.0 + prefactorA*q*q*q + prefactorB*q); };
                    splitting_function.setRange(0,1);
                    splitting_function.setTolerance(j.value("tab_utol",1e-9),j.value("tab_ftol",1e-2)); // Tolerance in splitting-function and its derivative
                    table = splitting_function.generate( PFF );
                }

                template<class Tparticle>
                    double operator()(const Tparticle &a, const Tparticle &b, double r2) {
                        if (r2 < rc2)
                            return Coulomb::operator()(a,b,r2)*splitting_function.eval(table,sqrt(r2)/rc1);
                        return 0;
                    }

                template<class Tparticle>
                    double operator()(const Tparticle &a, const Tparticle &b, const Point &r) const {
                        return operator()(a,b,r.squaredNorm());
                    }

                /**
                 * @param M2V Input of \f$ \frac{<M^2>}{9V\epsilon_0k_BT} \f$
                 * @brief Returns dielectric constant for the reaction field method (see DOI:10.1080/00268978300102721)
                 */
                double dielectric_constant(double M2V) const { 
                    if(eps_inf)
                        return (1.0 + 3.0*M2V);
                    if(eps_ins)
                        return ( 2.25*M2V + 0.25 + 0.75*sqrt(9.0*M2V*M2V + 2.0*M2V + 1.0) );
                    if(eps_vac)
                        return (2*M2V + 1.0)/(1.0 - M2V);
                    return 0.5*(2.0*eps_RF-1.0+sqrt(-72.0*M2V*M2V*eps_RF + 4.0*eps_RF*eps_RF + 4.0*eps_RF + 1.0))/(3*M2V-1.0); // Needs to be checked!
                    //return (6*M2V*eps_RF + 2*eps_RF + 1.0)/(1.0 + 2*eps_RF - 3*M2V); // Is OK when epsr=1.0
                }

                string info(char w) {
                    using namespace textio;
                    std::ostringstream o;
                    o << Coulomb::info(w)
                        << pad(SUB,w,"Cutoff") << rc1 << " "+angstrom << endl;
                    if(eps_inf) {
                        o << pad(SUB,w, epsilon_m+"_RF") << infinity << endl;
                    } else {
                        o << pad(SUB,w, epsilon_m+"_RF") << eps_RF << endl;
                    }
                    o << pad(SUB,w, "Inverse Debye length") << kappa << endl;
                    return o.str();
                }
        };

        /**
         * @brief Dipole-dipole interaction with reaction field.
         * 
         *  Keyword          |  Description
         * :--------------   | :---------------
         * `cutoff`          |  Cut-off for interactions.                             
         * `epsr`            |  Dielectric constant of the medium.                      (Default: 1)
         * `eps_rf`          |  Dielectric constant of the surroundings.            
         * 
         * @note If 'eps_rf' is set to (epsr,<0,0) then 'vacuum'/insulating/conducting boundary conditions will be used. 
         */
        class DipoleDipoleRF : public DipoleDipole {
            private:
                string _brief() { return "Dipole-dipole (reaction field)"; }
                double rc2,eps,eps_RF,epsr;
                bool eps_inf, eps_ins, eps_vac, eps_user; // Surrounding is: Conducting, Insulating, 'Vacuum', Set by user
            public:
                DipoleDipoleRF(Tmjson &j) : DipoleDipole(j) {
                    name+=" Reaction Field";
                    rc2 = pow(j.at("cutoff"),2.0);
                    epsr = j.value("epsr",1.0);
                    eps_RF = j.at("eps_rf");
                    eps_inf = false;
                    eps_ins = false;
                    eps_user = false;
                    eps_vac = false;
                    if(fabs(eps_RF) < 1e-6) {
                        eps_inf = true; // Conducting boundary conditions
                    } else if(eps_RF < 0.0) {
                        eps_ins = true; // Insulating boundary conditions
                    } else if(fabs(eps_RF-epsr) < 1e-6) {
                        eps_vac = true; // 'Vacuum' boundary conditions
                    } else {
                        eps_user = true; // Set by user
                    }
                    updateDiel(eps_RF);
                }

                template<class Tparticle>
                    double operator()(const Tparticle &a, const Tparticle &b, const Point &r) const {
                        if (r.squaredNorm() < rc2)
                            return (DipoleDipole::operator()(a,b,r) - eps*a.mu().dot(b.mu())*a.muscalar()*b.muscalar());
                        return 0.0;
                    }

                /** @brief Field at `r` due to dipole `p` 
                 * @warning Untested!
                 */
                template<class Tparticle>
                    Point field(const Tparticle &p, const Point &r) const {
                        return (DipoleDipole::field(p,r) + eps*p.mu()*p.muscalar());
                    }

                void updateDiel(double eps_rf_updated) {
                    if(eps_inf) {
                        eps = _lB/pow(rc2,1.5)/epsr;
                    } else if(eps_vac) {
                        eps = 0.0;
                    } else {
                        eps_RF = eps_rf_updated;
                        eps = _lB*(2*(eps_RF-epsr)/(2*eps_RF+epsr))/pow(rc2,1.5)/epsr;
                    }
                }  

                /**
                 * @param M2V Input of \f$ \frac{<M^2>}{9V\epsilon_0k_BT} \f$
                 * @brief Returns dielectric constant for the reaction field method (see DOI:10.1080/00268978300102721)
                 */
                double dielectric_constant(double M2V) const { 
                    if(eps_inf)
                        return (1.0 + 3.0*M2V);
                    if(eps_ins)
                        return ( 2.25*M2V + 0.25 + 0.75*sqrt(9.0*M2V*M2V + 2.0*M2V + 1.0) );
                    if(eps_vac)
                        return (2*M2V + 1.0)/(1.0 - M2V);
                    return 0.5*(2.0*eps_RF-1.0+sqrt(-72.0*M2V*M2V*eps_RF + 4.0*eps_RF*eps_RF + 4.0*eps_RF + 1.0))/(3*M2V-1.0); // Needs to be checked!
                    //return (6*M2V*eps_RF + 2*eps_RF + 1.0)/(1.0 + 2*eps_RF - 3*M2V); // Is OK when epsr=1.0
                }

                string info(char w) {
                    using namespace textio;
                    std::ostringstream o;
                    o << DipoleDipole::info(w)
                        << pad(SUB,w,"Cutoff") << sqrt(rc2) << " "+angstrom << endl;
                    if(eps_inf) {
                        o << pad(SUB,w, epsilon_m+"_RF") << infinity << endl;
                    } else {
                        o << pad(SUB,w, epsilon_m+"_RF") << eps_RF << endl;
                    }
                    return o.str();
                }
        };

        template<bool useIonIon=false, bool useIonDipole=false, bool useDipoleDipole=false, bool useIonQuadrupole=false>
            class MultipoleWolf : public PairPotentialBase {
                private:
                    WolfBase wolf;
                    string _brief() {
                        std::ostringstream o;
                        o << "Multipole Wolf, lB=" << _lB << textio::_angstrom;
                        return o.str();          
                    }
                protected:
                    double _lB;
                public:
                    //template<bool useIonIon=false, bool useIonDipole=false, bool useDipoleDipole=false, bool useIonQuadrupole=false>
                    MultipoleWolf(Tmjson &in) : wolf(in.value("kappa",0.0),in.at("cutoff")) {
                        name="Multipole Wolf";
                        _lB = Coulomb(in).bjerrumLength();
                    }
                    template<class Tparticle>
                        double operator()(const Tparticle &a, const Tparticle &b, const Point &r) {
                            double U_total = 0;
                            if((useIonIon? 1:0) + (useIonDipole? 1:0) + (useDipoleDipole? 1:0) + (useIonQuadrupole? 1:0) > 1) {
                                wolf.calcWolfData(r);
                                if (wolf.getR2i() < wolf.getRc2i())
                                    return 0;
                                if(useIonIon == true) U_total += wolf.q2q<true>(a.charge,b.charge,r);
                                if(useIonDipole == true) U_total += wolf.q2mu<true>(a.charge*b.muscalar(),b.mu(),b.charge*a.muscalar(),a.mu(),r);
                                if(useDipoleDipole == true) U_total += wolf.mu2mu<true>(a.mu(),b.mu(), a.muscalar()*b.muscalar(), r);
                                if(useIonQuadrupole == true) U_total += wolf.q2quad<true>(a.charge, b.theta,b.charge, a.theta,r);
                                return _lB*U_total;
                            }
                            if(useIonIon == true) U_total += wolf.q2q(a.charge,b.charge,r);
                            if(useIonDipole == true) U_total += wolf.q2mu(a.charge*b.muscalar(),b.mu(),b.charge*a.muscalar(),a.mu(),r);
                            if(useDipoleDipole == true) U_total += wolf.mu2mu(a.mu(),b.mu(), a.muscalar()*b.muscalar(), r);
                            if(useIonQuadrupole == true) U_total += wolf.q2quad(a.charge, b.theta,b.charge, a.theta,r);
                            return _lB*U_total;
                        }

                    template<bool useIon=true, bool useDipole=true, class Tparticle>
                        Point field(const Tparticle &p, const Point &r) {
                            if(useIon && useDipole) {
                                wolf.calcWolfData(r);
                                if (wolf.getR2i() < wolf.getRc2i())
                                    return Point(0,0,0);
                                Point E = wolf.fieldCharge<true>(p,r);
                                E += wolf.fieldDipole<true>(p,r);
                                return _lB*E;
                            }
                            if(useIon == true) return _lB*wolf.fieldCharge(p,r);
                            if(useDipole == true) return _lB*wolf.fieldDipole(p,r);
                            return Point(0,0,0);
                        }

                    /**
                     * @brief Interaction of dipole `p` with field `E`, see 'Intermolecular and SUrface Forces' by J. Israelachvili, p. 97 eq. 5.15
                     * @todo Needs to be tested!
                     */
                    template<class Tparticle>
                        double fieldEnergy(const Tparticle &p, const Point &E) {
                            return 0;
                        }

                    string info(char w) {
                        using namespace textio;
                        std::ostringstream o;
                        o << pad(SUB,w,"Temperature") << pc::T() << " K" << endl
                            << pad(SUB,w,"Bjerrum length") << _lB << " "+angstrom << endl
                            << pad(SUB,w,"Cutoff") << wolf.getCutoff() << " "+angstrom << endl
                            << pad(SUB,w,"Kappa") << wolf.getKappa() << " "+angstrom+"^-1" << endl;
                        return o.str();
                    }
            };

        class IonIonWolf : public Coulomb {
            private:
                string _brief() { return "Coulomb Wolf"; }
                WolfBase wolf;
            public:
                IonIonWolf(Tmjson &in) : Coulomb(in),
                wolf(in.at("kappa"), in.at("cutoff")) { 
                    name+=" Wolf"; 
                }

                template<class Tparticle>
                    double operator()(const Tparticle &a, const Tparticle &b, const Point &r) const {
                        return lB*wolf.q2q(a.charge,b.charge,r);
                    }

                template<class Tparticle>
                    Point field(const Tparticle &p, const Point &r) const {
                        return lB*wolf.fieldCharge(p,r);
                    }

                string info(char w) {
                    using namespace textio;
                    std::ostringstream o;
                    o << Coulomb::info(w)
                        << pad(SUB,w,"Cutoff") << wolf.getCutoff() << " "+angstrom << endl
                        << pad(SUB,w,"Kappa") << wolf.getKappa() << " "+angstrom+"^-1" << endl;
                    return o.str();
                }
        };

        class IonIonGaussianDamping : public Coulomb {
            private:
                string _brief() { return "Coulomb Gaussian Damping"; }
                GaussianDampingBase gdb;
            public:
                IonIonGaussianDamping(Tmjson &in) : Coulomb(in),
                gdb() { name+=" Gaussian Damping"; }

                template<class Tparticle>
                    double operator()(const Tparticle &a, const Tparticle &b, const Point &r) {
                        return lB*gdb.q2q(a.charge,b.charge,a.id,b.id,r);
                    }

                template<class Tparticle>
                    Point field(const Tparticle &p, const Point &r) const {
                        return lB*gdb.fieldCharge(p,r);
                    }
        };

        class IonDipoleGaussianDamping : public IonDipole {
            private:
                string _brief() { return "Ion-dipole Gaussian Damping"; }
                GaussianDampingBase gdb;
            public:
                IonDipoleGaussianDamping(Tmjson &in) : IonDipole(in),
                gdb() { name+=" Gaussian Damping"; }

                template<class Tparticle>
                    double operator()(const Tparticle &a, const Tparticle &b, const Point &r) {
                        return _lB*gdb.q2mu(a.charge*b.muscalar(),b.mu(),b.charge*a.muscalar(),a.mu(),a.id,b.id,r);
                    }
        };

        class DipoleDipoleGaussianDamping : public DipoleDipole {
            private:
                string _brief() { return "Dipole-dipole Gaussian Damping"; }
                GaussianDampingBase gdb;
            public:
                DipoleDipoleGaussianDamping(Tmjson &in) : DipoleDipole(in),
                gdb() { name+=" Gaussian Damping"; }

                template<class Tparticle>
                    double operator()(const Tparticle &a, const Tparticle &b, const Point &r) {
                        return _lB*gdb.mu2mu(a.mu(),b.mu(), a.muscalar()*b.muscalar(),a.id,b.id,r);
                    }

                template<class Tparticle>
                    Point field(const Tparticle &p, const Point &r) const {
                        return _lB*gdb.fieldDipole(p,r);
                    }
        };

        class IonQuadGaussianDamping : public IonQuad {
            private:
                string _brief() { return "Ion-Quadrupole Gaussian Damping"; }
                GaussianDampingBase gdb;
            public:
                IonQuadGaussianDamping(Tmjson &in) : IonQuad(in),
                gdb() { name+=" Gaussian Damping"; }

                template<class Tparticle>
                    double operator()(const Tparticle &a, const Tparticle &b, const Point &r) {
                        return _lB*gdb.q2quad(a.charge, b.theta,b.charge, a.theta,a.id,b.id,r);
                    }
        };

        /**
         * @brief Coulomb interaction with long-ranged compensation using radial-derivates 1-3 to be zero at the cut-off, see DOI: 10.1021/jp510612w.
         * The potential also mimics the interaction-tensor \f$ T_0(r,\alpha) = erfc(\alpha r)/r \f$ with \f$ \alpha = \sqrt{\pi} \f$.
         * 
         *  Keyword          |  Description
         * :--------------   | :---------------
         * `cutoff`          |  Cut-off for interactions.                               (Default: Infinity)
         * `eps_r`           |  Dielectric constant of the medium.                      (Default: \f$ \varepsilon_r = 1 \f$)
         * `tab_utol`        |  Tolerance in splitting-function error.                  (Default: \f$ 10^{-9}\f$)
         * `tab_ftol`        |  Tolerance of splitting-function derivative error.       (Default: \f$ 10^{-2}\f$)
         */
        class IonIonSP3 : public Coulomb {
            private:
                string _brief() { return "Coulomb SP3"; }
                double rc1, rc1i, rc2, _lB;
                Tabulate::Andrea<double> splitting_function;
                Tabulate::TabulatorBase<double>::data tabel;
            public:
                IonIonSP3(Tmjson &j) : Coulomb(j) { 
                    name += " SP3"; 
                    _lB = Coulomb(j).bjerrumLength();
                    rc1 = j.at("cutoff");
                    rc1i = 1.0/rc1;
                    rc2 = rc1*rc1;

                    std::function<double(double)> temp_function = [&](double q) { return (1.0 - 1.75*q + 5.25*q*q*q*q*q - 7.0*q*q*q*q*q*q + 2.5*q*q*q*q*q*q*q); };
                    splitting_function.setRange(0,rc1);
                    splitting_function.setTolerance(j.value("tab_utol",1e-9),j.value("tab_ftol",1e-2)); // Tolerance in splitting-function and its derivative
                    tabel = splitting_function.generate( temp_function );
                }

                template<class Tparticle>
                    double operator()(const Tparticle &a, const Tparticle &b, double r2) const {
                        double r1 = sqrt(r2);
                        if(r2 < rc2)
                            return _lB*(a.charge*b.charge/r1)*splitting_function.eval(tabel,r1*rc1i);
                        return 0.0;
                    }

                template<class Tparticle>
                    double operator()(const Tparticle &a, const Tparticle &b, const Point &r) const {
                        double r2 = r.squaredNorm();
                        return operator()(a,b,r2);
                    }

                /**
                 * @param M2V Input of \f$ \frac{<M^2>}{9V\epsilon_0k_BT} \f$
                 * @brief Returns dielectric constant for the SP3 method
                 */
                double dielectric_constant(double M2V) const { 
                    return (1.0 + 3.0*M2V);
                }

                string info(char w) {
                    using namespace textio;
                    std::ostringstream o;
                    o << Coulomb::info(w)
                        << pad(SUB,w,"Cutoff") << rc1 << " "+angstrom << endl
                        << splitting_function.info() << endl;
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

        /**
         * @brief Dipole-dipole interaction with long-ranged compensation using radial-derivates 1-3 to be zero at the cut-off for the original ionic interaction-tensor.
         * The potential also mimics the interaction-tensor \f$ T_0(r,\alpha) = erfc(\alpha r)/r \f$ with \f$ \alpha = \sqrt{\pi} \f$.
         * See DOI: 10.1021/jp510612w for more info.
         * 
         *  Keyword          |  Description
         * :--------------   | :---------------
         * `cutoff`          |  Cut-off for interactions.                               (Default: Infinity)
         * `eps_r`           |  Dielectric constant of the medium.                      (Default: \f$ \varepsilon_r = 1 \f$)
         */
        class DipoleDipoleSP3 : public DipoleDipole {
            private:
                string _brief() { return "DipoleDipole SP3"; }
                double rc1, rc1i, rc2, rc3;
            public:
                DipoleDipoleSP3(Tmjson &j) : DipoleDipole(j) { 
                    name += " SP3"; 
                    rc1 = j.at("cutoff");
                    rc1i = 1.0/rc1;
                    rc2 = rc1*rc1;
                    rc3 = rc2*rc1;
                }

                template<class Tparticle>
                    double operator()(const Tparticle &a, const Tparticle &b, const Point &r) const {
                        double r2 = r.squaredNorm();
                        if(r2 < rc2) {
                            double q = sqrt(r2)/rc1;
                            double q2 = q*q;
                            Eigen::Matrix3d T1 = r*r.transpose()/r2;
                            Eigen::Matrix3d T2 = Eigen::Matrix3d::Identity();
                            Eigen::Matrix3d T = ((42.0 - 105.0*q + 60.0*q2)*q2*T1 + (21.0 - 35.0*q + 15.0*q2)*q2*T2)/rc3;

                            double W = a.mu.transpose()*T*b.mu;
                            return (DipoleDipole::operator()(a,b,r) - _lB*W*a.muscalar*b.muscalar);
                        }
                        return 0.0;
                    }

                /**
                 * @param M2V Input of \f$ \frac{<M^2>}{9V\epsilon_0k_BT} \f$
                 * @brief Returns dielectric constant for the SP3 method
                 */
                double dielectric_constant(double M2V) const { 
                    return (1.0 + 3.0*M2V);
                }

                string info(char w) {
                    using namespace textio;
                    std::ostringstream o;
                    o << DipoleDipole::info(w)
                        << pad(SUB,w,"Cutoff") << rc1 << " "+angstrom << endl;
                    return o.str();
                }
        };

        /**
         * @brief Help-function for Qpotential
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
         * @brief Coulomb interaction with long-ranged compensation using moment cancellation, see Paper V in ISBN: 978-91-7422-440-5.
         *
         * Input parameters:
         * 
         *  Keyword    |  Description
         * :---------- | :---------------
         * `cutoff`    |  Cut-off distance for interactions (float, required)
         * `order`     |  Number of moments to cancel (default: 300)
         * `epsr`      |  Dielectric constant of the medium (default: \f$ \varepsilon_r = 1 \f$)
         * `tab_utol`  |  Tolerance in splitting-function error.                  (Default: \f$ 10^{-9}\f$)
         * `tab_ftol`  |  Tolerance of splitting-function derivative error.       (Default: \f$ 10^{-2}\f$)
         */
        template<class base=Coulomb, int k=1>
            class Qpotential : public base {
                private:
                    string _brief() { return base::brief() + ", Cutoff=" + std::to_string(rc1) + textio::_angstrom ; }
                    int order;
                    double rc1, rc1i, rc2;
                    Tabulate::Andrea<double> splitting_function;
                    Tabulate::TabulatorBase<double>::data tabel;
                public:
                    Qpotential(Tmjson &j) : base(j) { 
                        base::name += " Q"; 
                        rc1 = j.at("cutoff");
                        rc1i = 1.0/rc1;
                        rc2 = rc1*rc1;
                        order = j.value("order",300);

                        std::function<double(double)> Qk = [&](double q) { return qPochhammerSymbol(q,k,order); };  // Sets r^k-dependence, Ion-Ion -> k = 1, Ion-Dipole -> k = 2, etc. etc. 
                        splitting_function.setRange(0,1);
                        splitting_function.setTolerance(j.value("tab_utol",1e-9),j.value("tab_ftol",1e-2)); // Tolerance in splitting-function and its derivative
                        tabel = splitting_function.generate( Qk );
                    }

                    template<class Tparticle>
                        double operator()(const Tparticle &a, const Tparticle &b, double r2) const {
                            if ( r2 < rc2 )
                                return base::operator()(a, b, r2) * splitting_function.eval( tabel, sqrt(r2)  / rc1 );
                            return 0;
                        }

                    template<class Tparticle>
                        double operator()(const Tparticle &a, const Tparticle &b, const Point &r) const {
                            double r2 = r.squaredNorm();
                            if (r2 < rc2)
                                return base::operator()(a, b, r) * splitting_function.eval( tabel, sqrt(r2) / rc1 );
                            return 0;
                        }

                    /**
                     * @param M2V Input of \f$ \frac{<M^2>}{9V\epsilon_0k_BT} \f$
                     * @brief Returns dielectric constant for the ionic (k=1) or dipolar (k=3) \emph{q}-potential
                     */
                    double dielectric_constant(double M2V) const { 
                        if( k == 1 )
                            return (1.0 + 3.0*M2V);
                        if( k == 3 )
                            return (2*M2V + 1.0)/(1.0 - M2V);
                        return 0.0;
                    }

                    string info(char w) {
                        using namespace textio;
                        std::ostringstream o;
                        o << base::info(w)
                            << pad(SUB,w,"Cutoff") << rc1 << _angstrom << endl
                            << pad(SUB,w,"Order") << order << endl
                            << splitting_function.info() << endl;
                        return o.str();
                    }
            };

        /**
         * @brief Help-function for DipoleDipoleQ2.
         */
        inline double Euler_type_function(double q, int P=300, bool all=true) {
            if(q >= 1.0 - (1.0/2400.0))
                return 0.0;
            if(q <= (1.0/2400.0))
                return 1.0;
            double value1 = 0.0;
            double value2 = 0.0;
            double value3 = 0.0;
            for(int pt = -P; pt <= P; pt++) {
                double p = double(pt);
                if(all) {
                    value1 += pow(-1.0,p)*pow(q,(3.0*p*p-p)/2.0);
                    value2 += pow(-1.0,p)*(3.0*p*p-p)/(2.0*q)*pow(q,(3.0*p*p-p)/2.0);
                }
                value3 += pow(-1.0,p)*(3.0*p*p-p)/(2.0*q*q)*(3.0*p*p - p - 2.0)/2.0*pow(q,(3.0*p*p-p)/2.0);
            }
            return (value1 - value2*q + value3*q*q/3.0);
        }

        /**
         * @brief Dipole-dipole interaction with long-ranged compensation using moment cancellation.
         * This method is an expansion of the ion-ion interaction-tensor in Paper V in ISBN: 978-91-7422-440-5.
         * 
         *  Keyword          |  Description
         * :--------------   | :---------------
         * `cutoff`          |  Cut-off for interactions.                               (Default: Infinity)
         * `eps_r`           |  Dielectric constant of the medium.                      (Default: \f$ \varepsilon_r = 1 \f$)
         * `tab_utol`        |  Tolerance in prefactor-function error.                  (Default: \f$ 10^{-7}\f$)
         * `tab_ftol`        |  Tolerance of prefactor-function derivative error.       (Default: \f$ 10^{-2}\f$)
         * 
         * @warning Be careful with the splines, they tend to give diverging values for low `tab_utol` and `tab_ftol`!
         */
        class DipoleDipoleQ2 : public DipoleDipole {
            private:
                string _brief() { return "Dipole-dipole Q2"; }
                double rc1, rc1i, rc3i, tab_utol, tab_ftol;
                int N;
                Tabulate::Andrea<double> ak;
                Tabulate::Andrea<double> bk;
                Tabulate::TabulatorBase<double>::data tableA;
                Tabulate::TabulatorBase<double>::data tableB;
            public:
                DipoleDipoleQ2(Tmjson &j) : DipoleDipole(j) {
                    name += " Q2"; 
                    _lB = Coulomb(j).bjerrumLength();
                    rc1 = j.at("cutoff");
                    tab_utol = j.value("tab_utol",1e-7); // Higher accuracy than 1e-7 gives error
                    tab_ftol = j.value("tab_ftol",1e-2);
                    rc1i = 1.0/rc1;
                    rc3i = rc1i*rc1i*rc1i;
                    N = 300;

                    std::function<double(double)> Ak = [&](double q) { return Euler_type_function(q,N); };
                    ak.setRange(0,1);
                    ak.setTolerance(tab_utol,tab_ftol); // Tolerance in first prefactor-function and its derivative
                    tableA = ak.generate( Ak );

                    std::function<double(double)> Bk = [&](double q) { return Euler_type_function(q,N,false); };
                    bk.setRange(0.0,1.0);
                    bk.setTolerance(tab_utol,tab_ftol); // Tolerance in second prefactor-function and its derivative
                    tableB = bk.generate( Bk );
                }

                template<class Tparticle>
                    double operator()(const Tparticle &a, const Tparticle &b, const Point &r) const {
                        double r1 = r.norm();
                        if (r1 < rc1) {
                            double af = ak.eval(tableA,r1*rc1i);
                            double bf = bk.eval(tableB,r1*rc1i);
                            return _lB*mu2mu(a.mu, b.mu, a.muscalar*b.muscalar, r,af,bf);
                        }
                        return 0.0;
                    }

                /**
                 * @param M2V Input of \f$ \frac{<M^2>}{9V\epsilon_0k_BT} \f$
                 * @brief Returns dielectric constant for the expanded ionic \emph{q}-potential
                 */
                double dielectric_constant(double M2V) const { 
                    return (1.0 + 3.0*M2V);
                }

                string info(char w) {
                    using namespace textio;
                    std::ostringstream o;
                    o << DipoleDipole::info(w)
                        << pad(SUB,w,"Cutoff") << rc1 << " "+angstrom+"^-1" << endl;
                    o << ak.info() << endl;
                    o << bk.info() << endl;
                    return o.str();
                }
        };
    }
}
#endif

