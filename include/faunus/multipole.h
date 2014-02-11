#ifndef FAUNUS_MULTIPOLE_H
#define FAUNUS_MULTIPOLE_H
#include <faunus/common.h>
#include <faunus/auxiliary.h>
#include <faunus/species.h>
#include <faunus/picojson.h>

namespace Faunus {

  /**
   * @brief Approximation of erfc-function
   * @param x Value for which erfc should be calculated 
   * @details Reference for this approximation is found in Abramowitz and Stegun, Handbook of mathematical functions, eq. 7.1.26
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
   */
  double erfc_x(double x) {
    //
    // |error| <= 1.5*10^-7
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
  double erf_x(double x) {
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
   * @param QAxMuB Product of ion A:s charge and dipole B:s scalar
   * @param muB Unit dipole moment vector of particel B
   * @param QBxMuA Product of ion B:s charge and dipole A:s scalar
   * @param muA Unit dipole moment vector of particel A
   * @param r Direction \f$ r_{\mu_B} - r_{Q_A} \f$  // Check order!!!!!
   */
  template<class Tvec>
    double q2mu(double QAxMuB, const Tvec &muB, double QBxMuA, const Tvec &muA, const Tvec &r) {
      double r2i = 1/r.squaredNorm();
      double r1i = sqrt(r2i);
      double W1 = -QAxMuB*muB.dot(r);
      double W2 = QBxMuA*muA.dot(r);
      return (W1+W2)*r2i*r1i;  // Beware of r_Mu - r_Q = -r according to Israelachvili p.36, i.e. minus becomes plus
    }

  /**
   * @brief Returns dipole-dipole interaction
   *
   * @param muA Unit dipole moment vector of particle A
   * @param muB Unit dipole moment vector of particle B
   * @param muAxmuB Product of dipole scalars
   * @param r Vector \f$ r_{AB} \f$
   */
  template<class Tvec>
    double mu2mu(const Tvec &muA, const Tvec &muB, double muAxmuB, const Tvec &r) {
      double r2i = 1/r.squaredNorm();
      double r1i = sqrt(r2i);
      double r3i = r1i*r2i;
      double r5i = r3i*r2i;
      //Eigen::Matrix3d T = 3*r5i*r*r.transpose() - r3i*Matrix3d::Identity();
      //double W = -muA.transpose()*T*muB;                       // Buckingham    Å^-3
      double W = muA.dot(muB)*r3i-3*muA.dot(r)*muB.dot(r)*r5i; // J&K
      return W*muAxmuB;  // e^2 Å^2 Å ^-3 = e^2 /Å
    }

  /**
   * @brief Returns ion-quadrupole interaction
   * @param q Charge of particle A
   * @param quad Quadrupole moment of particle B
   * @param r Vector \f$ r_{AB} \f$
   */
  template<class Tvec, class Tmat>
    double q2quad(double q, const Tmat &quad, const Tvec &r) {
      double r2i = 1/r.squaredNorm();
      double r1i = sqrt(r2i);
      double r3i = r1i*r2i;
      double r5i = r3i*r2i;
      double W = r.transpose()*quad*r;
      W = W*r5i - quad.trace()*(r3i/3);
      return q*W; // e^2 / Å
    }




  /**
   * @brief Base class for Wolf based interactions
   *
   * The idea is that this class has no dependencies and is
   * to be used as a helper class for other classes.
   */
  class WolfBase {
    private:
      double rc1i_d, rc2i_d, rc3i_d, rc4i_d, rc5i_d, rc6i_d, rc1, rc1i, rc2i, rc3i, expKc, kappa, kappa2, kappa4, kappa6, constant;


      struct wdata {
        //template<class Tparticle>
        double r1i_d, r3i_d, r5i_d, diffR;
        //Tparticle* a,b;
      };

      template<class Tparticle, class Tvec>
        wdata calcWolfData(const Tparticle &a, const Tparticle &b, const Tvec &r) {
          wdata data;
          //data.a = a;
          //data.b = b;
          double r2i = 1/r.squaredNorm();
          double r1i = sqrt(r2i);
          double expK = 2*kappa*exp(-kappa2/r2i)/sqrt(pc::pi);
          data.diffR = (1/r1i) - (1/rc1i);
          data.r1i_d = erfc_x(kappa/r1i)*r1i;
          data.r3i_d = expK * (kappa2 + r2i) + data.r1i_d * r2i;
          data.r5i_d = expK*(r2i*r2i + (2*r2i/3)*kappa2 + (kappa6/(3*r2i)) + (kappa4/6)) + data.r1i_d*r2i*r2i;
          return data;
        }

    public:
      /**
       * @brief Constructor
       * @param alpha Dampening factor (inverse angstrom)
       * @param rcut Cutoff distance (angstrom)
       */
      WolfBase(double alpha, double rcut) {
        kappa = alpha;
        kappa2 = kappa*kappa;
        kappa4 = kappa2*kappa2;
        kappa6 = kappa4*kappa2;
        constant = 2/sqrt(pc::pi);
        rc1 = rcut;
        rc1i = 1/rcut;
        rc2i = rc1i*rc1i;
        rc3i = rc1i*rc2i;
        expKc = constant*kappa*exp(-kappa2/rc2i);
        rc1i_d = erfc_x(kappa/rc1i)*rc1i;
        rc2i_d = (expKc + rc1i_d)*rc1i;
        rc3i_d = expKc*(kappa2 + rc2i) + rc1i_d*rc2i;
        rc4i_d = expKc*( (2*kappa4/(3*rc1i)) + (2*kappa2*rc1i/3) + rc3i ) + rc1i_d*rc3i;
        rc5i_d = rc2i*expKc*(rc2i + (2/3)*kappa2 + (kappa6/(3*rc2i*rc2i)) + (kappa4/(6*rc2i))) + rc1i_d*rc2i*rc2i;
        rc6i_d = expKc*(rc2i*rc3i + (2/3)*kappa2*rc3i + (4/15)*kappa4*rc1i - (1/15)*(kappa6/rc1i) + (2/15)*(kappa4*kappa4/rc3i)) + rc1i_d*rc2i*rc3i;
      }

      /**
       * @brief Returns ion-dipole interaction, Needs to be checked!
       * @param QxMu Product of ion charge and dipole scalar
       * @param mu Unit dipole moment vector
       * @param r Direction \f$ r_Mu - r_Q \f$  
       */
      template<class Tvec>
        double q2mu(double QxMu1, const Tvec &mu1, double QxMu2, const Tvec &mu2, const Tvec &r) const {
          double r2i = 1/r.squaredNorm();
          if (r2i < rc2i)
            return 0;
          double r1i = sqrt(r2i);
          double expK = constant*kappa*exp(-kappa2/r2i);
          double r3i_d = expK * (kappa2 + r2i) + erfc_x(kappa/r1i)*r1i * r2i;
          double der = ((1/r1i) - rc1)*3*rc4i_d;
          double W1 = QxMu1*mu1.dot(r)*(r3i_d - rc3i_d + der);
          double W2 = QxMu2*mu2.dot(r)*(r3i_d - rc3i_d + der);
          return (W1-W2);  // Beware of r_Mu - r_Q = -r according to Israelachvili p.36, i.e. minus becomes plus
        }

      /**
       * @brief Dipole-dipole energy
       * @param muA Dipole moment (A) unit vector
       * @param muB Dipole moment (B) unit vector
       * @param muAxMuB Product of dipole moment scalars, |A|*|B|
       * @param r_ab Vector A to B
       * @param r1i inverse length of `r_ab`
       * @param r2i inverer length of `r_ab` squared
       * @returns energy in `kT/lB`
       */
      template<class Tvec>
        double mu2mu(const Tvec &muA, const Tvec &muB, double muAxmuB, const Tvec &r) const {
          double r1 = r.norm();
          double r2i = 1/r.squaredNorm();
          if (r2i < rc2i)
            return 0;
          double r1i = sqrt(r2i);
          double r1i_d = erfc_x(kappa/r1i)*r1i;
          double expK =  constant * kappa*exp(-kappa2/r2i);
          double r3i_d = expK * (kappa2 + r2i) + r1i_d * r2i;
          double r5i_d = r2i*expK*(r2i + (2/3)*kappa2 + (kappa6/(3*r2i*r2i)) + (kappa4/(6*r2i))) + r1i_d*r2i*r2i;
          double der = r1 - rc1;
          double t3 = muA.dot(muB)*(r3i_d - rc3i_d + der*3*rc4i_d);
          double t5 = 3*muA.dot(r)*muB.dot(r)*(r5i_d - rc5i_d + der*5*rc6i_d);
          //Eigen::Matrix3d T = Matrix3d::Identity()*(r3i_d - rc3i_d + der*3*rc4i_d) - 3*r*r.transpose()*(r5i_d - rc5i_d + der*5*rc6i_d);
          //double W = muA.transpose()*T*muB; 
          double W = t3 - t5;
          return W*muAxmuB;
        }

      /**
       * @brief Returns ion-quadrupole interaction
       */
      template<class Tvec, class Tmat>
        double q2quad(double q1, const Tmat &quad1, double q2, const Tmat &quad2, const Tvec &r) const {
          double r2i = 1/r.squaredNorm();
          if (r2i < rc2i)
            return 0;
          double r1i = sqrt(r2i);
          double expK =  constant * kappa*exp(-kappa2/r2i);
          double r1i_d = erfc_x(kappa/r1i)*r1i;
          double r3i_d = expK*(kappa2 + r2i) + r1i_d*r2i;
          double r5i_d = r2i*expK*(r2i + (2/3)*kappa2 + (kappa6/(3*r2i*r2i)) + (kappa4/(6*r2i))) + r1i_d*r2i*r2i;   
          double W1 = q1*r.transpose()*quad1*r;
          W1 = W1*(r5i_d - rc5i_d);
          W1 = - q1*quad1.trace()*((r3i_d - rc3i)/3);
          double W2 = q2*r.transpose()*quad2*r;
          W2 = W2*(r5i_d - rc5i_d);
          W2 = - q2*quad2.trace()*((r3i_d - rc3i_d)/3);
          return (W1+W2); // e^2 / Å
        }

      /** @brief Dipole field at `r` due to dipole `p` 
       *  Gets returned in [e/Å / lB] (\f$\beta eE / lB \f$)
       */
      template<class Tparticle>
        Point field(const Tparticle &p, const Tparticle &p0, const Point &r) const {
          double r2i = 1.0/r.squaredNorm();
          if (r2i < rc2i)
            return Point(0,0,0);
          double r1i = sqrt(r2i);
          Point r_n = r*r1i;
          double expK =  constant * kappa*exp(-kappa2/r2i);
          double r1i_d = erfc_x(kappa/r1i)*r1i;
          double der = (1/r1i) - rc1;
          //double r2i_d = expK*r1i + r1i_d*r1i;
          double r3i_d = expK*(kappa2 + r2i) + r1i_d*r2i;
          double r5i_d = expK*(r2i*r2i + (2/3)*kappa2*r2i + (kappa6/(3*r2i)) + (kappa4/6)) + r1i_d*r2i*r2i;
          Point fieldCharge = p.charge*r*r3i_d;
          Point fieldDipole = (3.0*p.mu.dot(r_n)*r_n*((r5i_d/r2i) - (rc5i_d/rc2i) + 3*der*(rc6i_d/rc2i)) - p.mu*(r3i_d - rc3i_d + 3*der*rc4i_d))*p.muscalar; // \beta e E
          return (fieldCharge + fieldDipole);
        }

      double getKappa() { return kappa; }
      double getCutoff() { return 1/rc1i; }
  };


























  class GaussianDamping {

    private:
      string _brief() { return "GaussianDamping"; }
    public:

      /**
       * @brief Returns ion-ion interaction, Needs to be checked!
       * @param QxQ Product of the ion charges
       * @param betaA Exponent of particles A:s Gaussian charge density
       * @param betaB Exponent of particles B:s Gaussian charge density
       * @param r Direction \f$ r_A - r_B \f$  
       */
      template<class Tvec>
        double q2q(double QxQ, double betaA, double betaB, const Tvec &r) {
          double beta12 = betaA*betaB/sqrt(betaA*betaA + betaB*betaB);
          double r1i = 1/(beta12*r.norm());
          return QxQ*beta12*erf_x(1/r1i)*r1i;
        }

      /**
       * @brief Returns ion-dipole interaction, Needs to be checked!
       * @param QAxMuB Product of ion charge and dipole scalar
       * @param muB Unit dipole moment vector
       * @param QBxMuA Product of ion charge and dipole scalar
       * @param muA Unit dipole moment vector
       * @param betaAC Exponent of particles A:s Gaussian charge density
       * @param betaBC Exponent of particles B:s Gaussian charge density
       * @param betaAD Exponent of particles A:s Gaussian dipole density
       * @param betaBC Exponent of particles B:s Gaussian dipole density
       * @param r Direction \f$ r_A - r_B \f$  
       */
      template<class Tvec>
        double q2mu(double QAxMuB, const Tvec &muB, double QBxMuA, const Tvec &muA, double betaAC, double betaBC, double betaAD, double betaBD, const Tvec &r) {
          double beta12 = betaAC*betaBD/sqrt(betaAC*betaAC + betaBD*betaBD);
          double r1i = 1/(beta12*r.norm());
          double B1A = r1i*r1i*(erf_x(1/r1i)*r1i - (2/sqrt(pc::pi))*exp(-1/(r1i*r1i)))*beta12*beta12*beta12;
          beta12 = betaAD*betaBC/sqrt(betaAD*betaAD + betaBC*betaBC);
          r1i = 1/(beta12*r.norm());
          double B1B = r1i*r1i*(erf_x(1/r1i)*r1i - (2/sqrt(pc::pi))*exp(-1/(r1i*r1i)))*beta12*beta12*beta12;
          double U1 = QAxMuB*muB.dot(r)*B1A;
          double U2 = QBxMuA*muA.dot(r)*B1B;
          return (U1 - U2);
        }

      /**
       * @brief Dipole-dipole energy, Needs to be checked!
       * @param muA Dipole moment (A) unit vector
       * @param muB Dipole moment (B) unit vector
       * @param muAxMuB Product of dipole moment scalars, |A|*|B|
       * @param betaA Exponent of particles A:s Gaussian charge density
       * @param betaB Exponent of particles B:s Gaussian charge density
       * @param r Direction \f$ r_A - r_B \f$ 
       * @returns energy in `kT/lB`
       */
      template<class Tvec>
        double mu2mu(const Tvec &muA, const Tvec &muB, double muAxmuB, double betaA, double betaB, const Tvec &r) const {
          double beta12 = betaA*betaB/sqrt(betaA*betaA + betaB*betaB);
          double beta123 = beta12*beta12*beta12;
          double r1i = 1/(beta12*r.norm()); 
          double r2i = r1i*r1i;
          double B1 = r2i*(erf_x(1/r1i)*r1i - (2/sqrt(pc::pi))*exp(-1/r2i));
          double B2 = 3*B1*r2i - (4/sqrt(pc::pi))*r2i*exp(-1/r2i);
          return muAxmuB*(muA.dot(muB)*beta123*B1 - (muA.dot(r))*(muB.dot(r))*beta123*beta123*B2);
        }

      /** 
       * @brief Field at `r` due to charge `p` 
       *  Gets returned in [e/Å / lB] (\f$\beta eE / lB \f$)
       * @param p Particle from where field is evaluated
       * @param beta Exponent of the Gaussian charge density of the particle where field is evaluated 
       * @param r Direction \f$ r_A - r_B \f$ 
       */
      template<class Tparticle>
        Point fieldQ2Q(const Tparticle &p,const Tparticle &p0, const Point &r) const {
          double beta12 = p.betaC*p0.betaC/sqrt(p.betaC*p.betaC + p0.betaC*p0.betaC);
          double r1i = 1/(beta12*r.norm()); 
          double r2i = r1i*r1i;
          double B1 = r2i*(erf_x(1/r1i)*r1i - (2/sqrt(pc::pi))*exp(-1/r2i));
          return (p.charge*beta12*beta12*beta12*B1*r);
        }

      /** 
       * @brief Field at `r` due to dipole `p` 
       *  Gets returned in [e/Å / lB] (\f$\beta eE / lB \f$)
       * @param p Particle from where field is evaluated
       * @param beta Exponent of the Gaussian dipole density of the particle where field is evaluated 
       * @param r Direction \f$ r_A - r_B \f$ 
       */
      template<class Tparticle>
        Point fieldMu2Mu(const Tparticle &p,const Tparticle &p0, const Point &r) const {
          double beta12 = p.betaD*p0.betaD/sqrt(p.betaD*p.betaD + p0.betaD*p0.betaD);
          double beta123 = beta12*beta12*beta12;
          double r1i = 1/(beta12*r.norm()); 
          double r2i = r1i*r1i;
          double B1 = r2i*(erf_x(1/r1i)*r1i - (2/sqrt(pc::pi))*exp(-1/r2i));
          double B2 = 3*B1*r2i - (4/sqrt(pc::pi))*r2i*exp(-1/r2i);
          return (beta123*B1*p.muscalar*p.mu - beta123*beta123*(p.muscalar*p.mu.dot(r))*B2*r);
        }

  };








  /**
   * @brief Base class for Wolf with Gaussian charge- or dipole-distribution based interactions
   *
   * The idea is that this class has no dependencies and is
   * to be used as a helper class for other classes.
   */

  class WolfGaussianDamping {
    private:
      double rc1, rc2i, constant, kappa, kappa2, kappa4, kappa6;
      Eigen::MatrixXd betaCC12, betaCC122, B0cCC, dB0cCC;
      Eigen::MatrixXd betaCD12, betaCD122, betaCD123, B1cCD, dB1cCD;
      Eigen::MatrixXd betaCQ12, betaCQ122, betaCQ123, B1cCQ, B2cCQ, dB1cCQ, dB2cCQ;
      Eigen::MatrixXd betaDD12, betaDD122, betaDD123, B1cDD, B2cDD, dB1cDD, dB2cDD;

      Eigen::VectorXd betaC, betaC2, betaC3, B1cC, dB1cC;
      Eigen::VectorXd betaD, betaD2, betaD3, B1cD, B2cD, dB1cD, dB2cD;  ////////////////////////////
      Eigen::VectorXd betaQ, betaQ2, betaQ3;
      
      struct wdata {
        //template<class Tparticle>
        double r2i, r1i_dW, r2i_dW, r3i_dW, r4i_dW, r5i_dW, erfX, expX, B0_s, B1_s, B2_s, diffR;
        //Tparticle* a,b;
      };

      template<class Tparticle, class Tvec>
        void calcWolfData(const Tparticle &a, const Tparticle &b, const Tvec &r, wdata &data) {
          data.r2i = 1/r.squaredNorm();
          double r1i = sqrt(data.r2i);
          double r3i = data.r2i*r1i;
          double r4i = data.r2i*data.r2i;
          double expK =  constant*kappa*exp(-kappa2/data.r2i);
          data.r1i_dW = erfc_x(kappa/r1i)*r1i;
          data.r2i_dW = (expK + data.r1i_dW)*r1i;
          data.r3i_dW = expK*(kappa2 + data.r2i) + data.r1i_dW * data.r2i;
          data.r4i_dW = expK*((2*kappa4/(3*r1i)) + (2*kappa2*r1i/3) + r3i) + data.r1i_dW*r3i;
          data.r5i_dW = expK*data.r2i*(data.r2i + (2/3)*kappa2 + (kappa6/(3*r4i)) + (kappa4/(6*data.r2i))) + data.r1i_dW*r4i;   
          data.erfX = erf_x(betaDD12(a.id-1,b.id-1)/r1i);
          data.expX = constant*betaDD12(a.id-1,b.id-1)*exp(-betaDD12(a.id-1,b.id-1)/data.r2i);

          data.B0_s = erf_x(betaCC12(a.id-1,b.id-1)/r1i)*data.r1i_dW;
          data.B1_s = data.erfX*data.r3i_dW - data.expX*data.r2i_dW;            
          data.B2_s = 3*data.erfX*data.r5i_dW - data.expX*(2*data.r2i_dW*betaDD122(a.id-1,b.id-1) + 3*data.r4i_dW);
          data.diffR = (1/r1i) - rc1;
          return data;
        }

    public:
      /**
       * @brief Constructor
       * @param p Particles
       * @param kappa_in Dampening factor (inverse angstrom)
       * @param betaC Exponent of Gaussian charge distribution
       * @param betaD Exponent of Gaussian dipole distribution
       * @param rcut Cutoff distance (angstrom)
       */
      WolfGaussianDamping(double kappa_in, double rcut) {
        //WolfGaussianDamping(const Tparticles &p, double kappa_in, double rcut) {
        Eigen::MatrixXd expBCC, expBCD, expBDD, expBCQ, erfBCC, erfBCD, erfBDD, erfBCQ;
        Eigen::VectorXd expBC, expBD, expBQ, erfBC, erfBD, erfBQ;
        constant = 2/sqrt(pc::pi);
        kappa = kappa_in;
        kappa2 = kappa*kappa;
        kappa4 = kappa2*kappa2;
        kappa6 = kappa4*kappa2;
        rc1 = rcut;
        double rc1i = 1/rc1;
        rc2i = rc1i*rc1i;
        double rc3i = rc2i*rc1i;
        double expKc = constant*kappa*exp(-kappa2/rc2i);

        double rc1i_dW = erfc_x(kappa/rc1i)*rc1i;
        double rc2i_dW = (expKc + rc1i_dW)*rc1i;
        double rc3i_dW = expKc*(kappa2 + rc2i) + rc1i_dW*rc2i;
        double rc4i_dW = expKc*((2*kappa4/(3*rc1i)) + (2*kappa2*rc1i/3) + rc3i) + rc1i_dW*rc3i;
        double rc5i_dW = expKc*(rc2i*rc2i + (2*rc2i*kappa2/3) + (kappa6/(3*rc2i)) + (kappa4/6)) + rc1i_dW*rc2i*rc2i;
        double rc6i_dW = expKc*(rc2i*rc3i + (2/3)*kappa2*rc3i + (4/15)*kappa4*rc1i - (1/15)*(kappa6/rc1i) + (2/15)*(kappa4*kappa4/rc3i)) + rc1i_dW*rc2i*rc3i;

        int N = atom.size() - 1;
        betaC.resize(N);
        betaD.resize(N);
        betaQ.resize(N);
        betaC2.resize(N);
        betaD2.resize(N);
        betaQ2.resize(N);
        betaC3.resize(N);
        betaD3.resize(N);
        betaQ3.resize(N);
        erfBC.resize(N);
        erfBD.resize(N);
        erfBQ.resize(N);
        expBC.resize(N);
        expBD.resize(N);
        expBQ.resize(N);
        B1cC.resize(N);
        B1cD.resize(N);
        B2cD.resize(N);
        dB1cC.resize(N);
        dB1cD.resize(N);
        dB2cD.resize(N);
        erfBCC.resize(N,N);
        erfBCD.resize(N,N);
        erfBDD.resize(N,N);
        erfBCQ.resize(N,N);
        expBCC.resize(N,N);
        expBCD.resize(N,N);
        expBDD.resize(N,N);
        expBCQ.resize(N,N);
        betaCC12.resize(N,N);
        betaCD12.resize(N,N);
        betaCQ12.resize(N,N);
        betaCC122.resize(N,N);
        betaCD122.resize(N,N);
        betaCD123.resize(N,N);
        betaDD12.resize(N,N);
        betaDD122.resize(N,N);
        betaCQ122.resize(N,N);
        betaDD123.resize(N,N);
        betaCQ123.resize(N,N);
        B0cCC.resize(N,N);
        B1cCD.resize(N,N);
        B1cDD.resize(N,N);
        B1cCQ.resize(N,N);
        B2cDD.resize(N,N);
        B2cCQ.resize(N,N);
        dB0cCC.resize(N,N);
        dB1cCD.resize(N,N);
        dB1cDD.resize(N,N);
        dB1cCQ.resize(N,N);
        dB2cDD.resize(N,N);
        dB2cCQ.resize(N,N);

        for(int i = 0; i < N; i++) {
          betaC(i) = atom[i+1].betaC;
          betaC2(i) = betaC(i)*betaC(i);
          betaC3(i) = betaC2(i)*betaC(i);
          betaD(i) = atom[i+1].betaD;
          betaD2(i) = betaD(i)*betaD(i);
          betaD3(i) = betaD2(i)*betaD(i);
          betaQ(i) = atom[i+1].betaQ;
          betaQ2(i) = betaQ(i)*betaQ(i);
          betaQ3(i) = betaQ2(i)*betaQ(i);
        }

        for(int i = 0; i < N; i++) {
          for(int j = 0; j < N; j++) {
            betaCC12(i,j) = betaC(i)*betaC(j)/sqrt(betaC(i)*betaC(i) + betaC(j)*betaC(j));
            betaCD12(i,j) = betaC(i)*betaD(j)/sqrt(betaC(i)*betaC(i) + betaD(j)*betaD(j));
            betaDD12(i,j) = betaD(i)*betaD(j)/sqrt(betaD(i)*betaD(i) + betaD(j)*betaD(j));
            betaCQ12(i,j) = betaC(i)*betaQ(j)/sqrt(betaC(i)*betaC(i) + betaQ(j)*betaQ(j));
            betaCC122(i,j) = betaCC12(i,j)*betaCC12(i,j);
            betaCD122(i,j) = betaCD12(i,j)*betaCD12(i,j);
            betaDD122(i,j) = betaDD12(i,j)*betaDD12(i,j);
            betaCQ122(i,j) = betaCQ12(i,j)*betaCQ12(i,j);
            betaCD123(i,j) = betaCD122(i,j)*betaCD12(i,j);
            betaDD123(i,j) = betaDD122(i,j)*betaDD12(i,j);
            betaCQ123(i,j) = betaCQ122(i,j)*betaCQ12(i,j);

            expBCC(i,j) = constant*betaCC12(i,j)*exp(-betaCC122(i,j)/rc2i);
            expBCD(i,j) = constant*betaCD12(i,j)*exp(-betaCD122(i,j)/rc2i);
            expBDD(i,j) = constant*betaDD12(i,j)*exp(-betaDD122(i,j)/rc2i);
            expBCQ(i,j) = constant*betaCQ12(i,j)*exp(-betaCQ122(i,j)/rc2i);
            erfBCC(i,j) = erf_x(betaCC12(i,j)/rc1i);
            erfBCD(i,j) = erf_x(betaCD12(i,j)/rc1i);
            erfBDD(i,j) = erf_x(betaDD12(i,j)/rc1i);
            erfBCQ(i,j) = erf_x(betaCQ12(i,j)/rc1i);

            B0cCC(i,j) = erfBCC(i,j)*rc1i_dW;
            B1cCD(i,j) = erfBCD(i,j)*rc3i_dW - expBCD(i,j)*rc2i_dW;
            B1cDD(i,j) = erfBDD(i,j)*rc3i_dW - expBDD(i,j)*rc2i_dW;    
            B1cCQ(i,j) = erfBCQ(i,j)*rc3i_dW - expBCQ(i,j)*rc2i_dW; 
            B2cDD(i,j) = 3*erfBDD(i,j)*rc5i_dW - expBDD(i,j)*(2*rc2i_dW*betaDD122(i,j) + 3*rc4i_dW);
            B2cCQ(i,j) = 3*erfBCQ(i,j)*rc5i_dW - expBCQ(i,j)*(2*rc2i_dW*betaCQ122(i,j) + 3*rc4i_dW);
            
            dB0cCC(i,j) = expBCC(i,j)*rc1i_dW - erfBCC(i,j)*rc2i_dW;
            dB1cCD(i,j) = expBCD(i,j)*(rc3i_dW + (2/betaCD12(i,j))*rc3i_dW + 2*betaCD12(i,j)*rc1i_dW) - 3*erfBCD(i,j)*rc4i_dW;
            dB1cDD(i,j) = expBDD(i,j)*(rc3i_dW + (2/betaDD12(i,j))*rc3i_dW + 2*betaDD12(i,j)*rc1i_dW) - 3*erfBDD(i,j)*rc4i_dW;
            dB1cCQ(i,j) = expBCQ(i,j)*(rc3i_dW + (2/betaCQ12(i,j))*rc3i_dW + 2*betaCQ12(i,j)*rc1i_dW) - 3*erfBCQ(i,j)*rc4i_dW;
            dB2cDD(i,j) = expBDD(i,j)*(15*rc5i_dW + 10*betaDD122(i,j)*rc3i_dW + 4*betaDD122(i,j)*betaDD122(i,j)*rc1i_dW) - 15*erfBDD(i,j)*rc6i_dW;
            dB2cCQ(i,j) = expBCQ(i,j)*(15*rc5i_dW + 10*betaCQ122(i,j)*rc3i_dW + 4*betaCQ122(i,j)*betaCQ122(i,j)*rc1i_dW) - 15*erfBCQ(i,j)*rc6i_dW;
          }
          
          expBC(i) = constant*betaC(i)*exp(-betaC2(i)/rc2i);
          expBD(i) = constant*betaD(i)*exp(-betaD2(i)/rc2i);
          expBQ(i) = constant*betaQ(i)*exp(-betaQ2(i)/rc2i);
          erfBC(i) = erf_x(betaC(i)/rc1i);
          erfBD(i) = erf_x(betaD(i)/rc1i);
          erfBQ(i) = erf_x(betaQ(i)/rc1i);
          
          B1cC(i) = erfBC(i)*rc3i_dW - expBC(i)*rc2i_dW;
          B1cD(i) = erfBD(i)*rc3i_dW - expBD(i)*rc2i_dW;    
          B2cD(i) = 3*erfBD(i)*rc5i_dW - expBD(i)*(2*rc2i_dW*betaD2(i) + 3*rc4i_dW);
            
          dB1cC(i) = expBC(i)*(rc3i_dW + (2/betaC(i))*rc3i_dW + 2*betaC(i)*rc1i_dW) - 3*erfBC(i)*rc4i_dW;
          dB1cD(i) = expBD(i)*(rc3i_dW + (2/betaD(i))*rc3i_dW + 2*betaD(i)*rc1i_dW) - 3*erfBD(i)*rc4i_dW;
          dB2cD(i) = expBD(i)*(15*rc5i_dW + 10*betaD2(i)*rc3i_dW + 4*betaD2(i)*betaD2(i)*rc1i_dW) - 15*erfBD(i)*rc6i_dW;
        }
      }

      /**
       * @brief Returns ion-ion interaction. Needs to be checked!
       * @param QxQ Product of ion charges
       * @param betaA Index of type atom A:s charge distribution
       * @param betaB Index of type atom B:s charge distribution
       * @param r Direction \f$ r_{q_A} - r_{q_B} \f$  
       */
      template<class Tvec>
        double q2q(double QxQ, int betaA, int betaB, const Tvec &r) {
          double r1 = r.norm();
          if(r1 > rc1)
            return 0;
          double r1i_dW = erfc_x(kappa*r1)/r1;
          double B0 = erf_x(betaCC12(betaA,betaB)*r1)*r1i_dW;
          return QxQ*(B0 - B0cCC(betaA,betaB) - (r1 - rc1)*dB0cCC(betaA,betaB));
        }

      /**
       * @brief Returns ion-dipole interaction, Needs to be checked!
       * @param QBxMuA Product of ion B:s charge and dipole A:s scalar
       * @param muA Unit dipole moment vector of particle A
       * @param QAxMuB Product of ion A:s charge and dipole B:s scalar
       * @param muB Unit dipole moment vector of particle B
       * @param betaA Index of type atom A:s charge/dipole distribution
       * @param betaB Index of type atom B:s dipole/charge distribution
       * @param r Direction \f$ r_A - r_B \f$  
       */
      template<class Tvec>
        double q2mu(double QBxMuA, const Tvec &muA, double QAxMuB, const Tvec &muB, int betaA, int betaB, const Tvec &r) const {
          double r2i = 1/r.squaredNorm();
          if (r2i < rc2i)
            return 0;
          double r1i = sqrt(r2i);
          double expK = constant*kappa*exp(-kappa2/r2i);
          double r1i_dW = erfc_x(kappa/r1i)*r1i;
          double r2i_dW = (expK + r1i_dW)*r1i;
          double r3i_dW = expK*(kappa2 + r2i) + r1i_dW*r2i;
          double B1 = erf_x(betaCD12(betaA,betaB)/r1i)*r3i_dW - constant*betaCD12(betaA,betaB)*exp(-betaCD122(betaA,betaB)/r2i)*r2i_dW;
          return (QBxMuA*muA.dot(r) - QAxMuB*muB.dot(r))*(B1 - B1cCD(betaA,betaB) - ((1/r1i) - rc1)*dB1cCD(betaA,betaB));  // Beware of r_Mu - r_Q = -r according to Israelachvili p.36, i.e. minus becomes plus
        }

      /**
       * @brief Dipole-dipole energy
       * @param muA Dipole moment (A) unit vector
       * @param muB Dipole moment (B) unit vector
       * @param muAxMuB Product of dipole moment scalars, |A|*|B|
       * @param betaA Index of type atom A:s dipole distribution
       * @param betaB Index of type atom B:s dipole distribution
       * @param r_ab Vector A to B
       * @returns energy in `kT/lB`
       */
      template<class Tvec>
        double mu2mu(const Tvec &muA, const Tvec &muB, double muAxmuB, int betaA, int betaB, const Tvec &r) const {
          double r2i = 1/r.squaredNorm();
          if (r2i < rc2i)
            return 0;
          double r1i = sqrt(r2i);
          double der = (1/r1i) - rc1;
          double r3i = r2i*r1i;
          double expK =  constant*kappa*exp(-kappa2/r2i);
          double r1i_dW = erfc_x(kappa/r1i)*r1i;
          double r2i_dW = (expK + r1i_dW)*r1i;
          double r3i_dW = expK*(kappa2 + r2i) + r1i_dW * r2i;
          double r4i_dW = expK*((2*kappa4/(3*r1i)) + (2*kappa2*r1i/3) + r3i) + r1i_dW*r3i;
          double r5i_dW = expK*(r2i*r2i + (2*r2i*kappa2/3) + (kappa6/(3*r2i)) + (kappa4/6)) + r1i_dW*r2i*r2i;  
          double erfX = erf_x(betaDD12(betaA,betaB)/r1i);
          double expX = constant*betaDD12(betaA,betaB)*exp(-betaDD122(betaA,betaB)/r2i);

          double B1_s = erfX*r3i_dW - expX*r2i_dW;            
          double B2_s = 3*erfX*r5i_dW - expX*(2*r2i_dW*betaDD122(betaA,betaB) + 3*r4i_dW);
          double W = muA.dot(muB)*(B1_s - B1cDD(betaA,betaB) - der*dB1cDD(betaA,betaB)) - muA.dot(r)*muB.dot(r)*(B2_s - B2cDD(betaA,betaB) -der*dB2cDD(betaA,betaB));
          return W*muAxmuB;
        }

        /**
         * @brief Returns ion-quadrupole interaction. Not implemented!
         */
        
        template<class Tvec, class Tmat>
          double q2quad(double q1, const Tmat &quad1, double q2, const Tmat &quad2, int betaA, int betaB, const Tvec &r) const {
            double r2i = 1/r.squaredNorm();
            if (r2i < rc2i)
              return 0;
            double r1i = sqrt(r2i);
            double der = (1/r1i) - rc1;
            double r3i = r1i*r2i;
            double expK =  2 * kappa*exp(-kappa2/r2i) / sqrt(pc::pi);
            double r1i_dW = erfc_x(kappa/r1i)*r1i;
            double r2i_dW = (expK + r1i_dW)*r1i;
            double r3i_dW = expK*(kappa2 + r2i) + r1i_dW*r2i;
            double r4i_dW = expK*((2*kappa4/(3*r1i)) + (2*kappa2*r1i/3) + r3i) + r1i_dW*r3i;
            double r5i_dW = r2i*expK*(r2i + (2/3)*kappa2 + (kappa6/(3*r2i*r2i)) + (kappa4/(6*r2i))) + r1i_dW*r2i*r2i;  
            double erfX = erf_x(betaDD12(betaA,betaB)/r1i);
            double expX = constant*betaDD12(betaA,betaB)*exp(-betaDD122(betaA,betaB)/r2i);

            double B1_s = erfX*r3i_dW - expX*r2i_dW;            
            double B2_s = 3*erfX*r5i_dW - expX*(2*r2i_dW*betaDD122(betaA,betaB) + 3*r4i_dW);

            Eigen::MatrixXd Temp = q1*quad1 + q2*quad2;
            double W = r.transpose()*Temp*r;
            W = W*(B2_s - B2cCQ(betaA,betaB) - der*dB2cCQ(betaA,betaB));
            W = W  - (q1*quad1.trace() + q2*quad2.trace())*(B1_s - B1cCQ(betaA,betaB) -der*dB1cCQ(betaA,betaB));
            return (W/3); // e^2 / Å
          }

        /** 
         * @brief Dipole field at `r` due to dipole `p` 
         * @param p Particle from which to evaluate the field
         * @param betaA Index of type atom A:s dipole distribution
         * @param betaB Index of type atom B:s dipole distribution
         * @param r_ab Vector A to B
         * @returns [e/Å / lB] (\f$\beta eE / lB \f$)
         */
        template<class Tparticle>
          Point field(const Tparticle &p, int beta, const Point &r) const {
            double r2i = 1/r.squaredNorm();
            if (r2i < rc2i)
              return Point(0,0,0);
            double r1i = sqrt(r2i);
            double r3i = r2i*r1i;
            double der = (1/r1i) - rc1;
            double expK = constant*kappa*exp(-kappa2/r2i);
            double r1i_dW = erfc_x(kappa/r1i)*r1i;
            double r2i_dW = (expK + r1i_dW)*r1i;
            double r3i_dW = expK*(kappa2 + r2i) + r1i_dW*r2i;
            double r4i_dW = expK*((2*kappa4/(3*r1i)) + (2*kappa2*r1i/3) + r3i) + r1i_dW*r3i;
            double r5i_dW = expK*r2i*(r2i + (2/3)*kappa2 + (kappa6/(3*r2i*r2i)) + (kappa4/(6*r2i))) + r1i_dW*r2i*r2i;
            double erfX = erf_x(betaD(beta)/r1i);
            double expX = constant*exp(-betaD2(beta)/r2i);
            double B1_c = erf_x(betaC(beta)/r1i)*r3i_dW - constant*betaC3(beta)*exp(-betaC2(beta)/r2i)*r2i_dW;
            double B1_d = erf_x(betaD(beta)/r1i)*r3i_dW - constant*betaD3(beta)*exp(-betaD2(beta)/r2i)*r2i_dW;
            double B2 = 3*erfX*betaD(beta)*r5i_dW - expX*betaD2(beta)*(2*r2i_dW*betaD2(beta) + 3*r4i_dW);
            Point fieldCharge = p.charge*r*(B1_c - B1cC(beta) - ((1/r1i) - rc1)*dB1cC(beta));
            Point fieldDipole = -p.muscalar*(p.mu*(B1_d - B1cD(beta) - der*dB1cD(beta)) - r*p.mu.dot(r)*(B2 - B2cD(beta) -der*dB2cD(beta)));
            return (fieldCharge + fieldDipole);
          }

          double getKappa() { return kappa; }
          double getCutoff() { return rc1; }
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
              IonDipole(InputMap &in) {
                pc::setT ( in.get<double>("temperature", 298.15, "Absolute temperature (K)") );
                double epsilon_r = in.get<double>("epsilon_r",80., "Dielectric constant");
                _lB=pc::lB( epsilon_r );
              }
              /**
               * @brief Ion-dipole
               * @param a Dipole particle A
               * @param b Dipole particle B
               * @param r Direction \f$ r_A - r_B \f$  
               */
              template<class Tparticle> // q2mu(1->2,r) + q2mu(2->1,-r) = q2mu(1->2,r) - q2mu(2->1,r)
                double operator()(const Tparticle &a, const Tparticle &b, const Point &r) const {
                  return _lB*q2mu(a.charge*b.muscalar,b.mu,b.charge*a.muscalar,a.mu,r);
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
              DipoleDipole(InputMap &in) {
                name="Dipole-dipole";
                pc::setT ( in.get<double>("temperature", 298.15,
                      "Absolute temperature (K)") );
                double epsilon_r = in.get<double>("epsilon_r",80.,
                    "Dielectric constant");
                _lB = pc::lB(epsilon_r);
              }
              template<class Tparticle>
                double operator()(const Tparticle &a, const Tparticle &b, const Point &r) const {
                  return _lB*mu2mu(a.mu, b.mu, a.muscalar*b.muscalar, r);
                }

              /** @brief Dipole field at `r` due to dipole `p` 
               *  Gets returned in [e/Å] (\f$\beta eE \f$)
               */
              template<class Tparticle>
                Point field(const Tparticle &p,const Tparticle &p0, const Point &r) const {
                  double R2 = 1.0/r.squaredNorm();
                  double R1 = sqrt(R2);
                  Point r_n = r*R1;
                  return ((3.0*p.mu.dot(r_n)*r_n - p.mu)*R2*R1)*p.muscalar*_lB; // \beta e E
                }

              /**
               * @brief Interaction of dipole `p` with field `E`, see 'Intermolecular and SUrface Forces' by J. Israelachvili, p. 97 eq. 5.15
               * @todo Needs to be tested!
               */
              template<class Tparticle>
                double fieldEnergy(const Tparticle &p, const Point &E) {
                  return -p.muscalar*p.mu.dot(E);
                }

              string info(char w) {
                using namespace textio;
                std::ostringstream o;
                o << pad(SUB,w,"Temperature") << pc::T() << " K" << endl
                  << pad(SUB,w,"Bjerrum length") << _lB << " "+angstrom << endl;
                return o.str();
              }
          };

          /**
           * @brief Dipole-dipole interaction w. spherical cutoff and reaction field
           *
           * More info...
           */
          class DipoleDipoleRF : public DipoleDipole {
            private:
              string _brief() { return "Dipole-dipole (RF)"; }
              double rc2,eps,eps_rf;
            public:
              DipoleDipoleRF(InputMap &in) : DipoleDipole(in) {
                name+=" Reaction Field";
                rc2 = pow(in.get<double>("dipdip_cutoff",pc::infty), 2);
                eps_rf = in.get<double>("epsilon_rf",80.);
                eps = _lB*(2*(eps_rf-1)/(2*eps_rf+1))/pow(rc2,1.5);
              }
              template<class Tparticle>
                double operator()(const Tparticle &a, const Tparticle &b, const Point &r) const {
                  if (r.squaredNorm() < rc2)
                    return (DipoleDipole::operator()(a,b,r) - eps*a.mu.dot(b.mu)*a.muscalar*b.muscalar);
                  return 0;
                }

              void updateDiel(double er) {
                eps = _lB*(2*(er-1)/(er+1))/pow(rc2,1.5);
              }  

              string info(char w) {
                using namespace textio;
                std::ostringstream o;
                o << DipoleDipole::info(w)
                  << pad(SUB,w,"Cutoff") << sqrt(rc2) << " "+angstrom << endl
                  << pad(SUB,w,"epsilon_rf") << eps_rf << endl;
                return o.str();
              }
          };

          /**
           * @brief Ion-dipole interaction
           *
           * More info...
           */
          class IonQuad : public PairPotentialBase {
            private:
              string _brief() { return "Ion-quadrupole"; }
            protected:
              double _lB;
            public:
              IonQuad(InputMap &in) {
                name="Ion-Quad";
                pc::setT ( in.get<double>("temperature", 298.15, "Absolute temperature (K)") );
                double epsilon_r = in.get<double>("epsilon_r",80., "Dielectric constant");
                _lB=pc::lB( epsilon_r );
              }
              template<class Tparticle>
                double operator()(const Tparticle &a, const Tparticle &b, const Point &r) const {
                  return _lB*(q2quad(a.charge, b.theta,r)+q2quad(b.charge, a.theta,r));
                }
                
              /** @brief Quadrupole field at `r` due to quadrupole `p` 
               *  Gets returned in [e/Å] (\f$\beta eE \f$)
               */
              template<class Tparticle>
                Point field(const Tparticle &p,const Tparticle &p0, const Point &r) const {
                  double r2i = 1.0/r.squaredNorm();
                  double r1i = sqrt(r2i);
                  double r4i = r2i*r2i;
                  Point temp1 = r.transpose()*p.theta;
                  double temp2 = temp1.dot(r);
                  temp1 = p.theta*r;
                  Point E = (temp2*5*r4i*r2i - p.theta.trace()*r4i)*r - 2*temp1*r4i*r1i; // Check 1/r dep
                  return E;
                }

              string info(char w) { return _brief(); }
          };

          class IonDipoleWolf : public IonDipole {
            private:
              string _brief() { return "Ion-dipole Wolf"; }
              WolfBase wolf;
            public:
              IonDipoleWolf(InputMap &in) : IonDipole(in),
              wolf(in.get<double>("kappa", 1.8, "Kappa-damping"),
                  in.get<double>("dipdip_cutoff",in.get<double>("cuboid_len",pc::infty)/2)) {
                name+=" Wolf";
              }

              template<class Tparticle>
                double operator()(const Tparticle &a, const Tparticle &b, const Point &r) const {
                  return _lB*wolf.q2mu(b.charge*a.muscalar, a.mu, a.charge*b.muscalar, b.mu,r);
                }

              template<class Tparticle>
                Point field(const Tparticle &p,const Tparticle &p0, const Point &r) const {
                  return Point(0,0,0);
                }

              string info(char w) {
                using namespace textio;
                std::ostringstream o;
                o << IonDipole::info(w)
                  << pad(SUB,w,"Cutoff") << wolf.getCutoff() << " "+angstrom << endl
                  << pad(SUB,w,"Kappa") << wolf.getKappa() << " "+angstrom+"^-1" << endl;
                return o.str();
              }
          };
          
          class IonDipoleWolfDamped : public IonDipole {
            private:
              string _brief() { return "Ion-dipole Wolf Damped"; }
              WolfGaussianDamping wolf;
            public:
              IonDipoleWolfDamped(InputMap &in) : IonDipole(in),
              wolf(in.get<double>("kappa", 1.8, "Kappa-damping"),
                  in.get<double>("dipdip_cutoff",in.get<double>("cuboid_len",pc::infty)/2)) {
                name+=" Wolf Damped";
              }

              template<class Tparticle>
                double operator()(const Tparticle &a, const Tparticle &b, const Point &r) const {
                  return _lB*wolf.q2mu(b.charge*a.muscalar, a.mu, a.charge*b.muscalar, b.mu, a.id-1, b.id-1,r);
                }

              /** @brief Dipole field at `r` due to dipole `p` 
               *  Gets returned in [e/Å] (\f$\beta eE \f$)
               */
              template<class Tparticle>
                Point field(const Tparticle &p,const Tparticle &p0, const Point &r) const {
                  return _lB*wolf.field(p,p0,r);
                }

              string info(char w) {
                using namespace textio;
                std::ostringstream o;
                o << IonDipole::info(w)
                  << pad(SUB,w,"Cutoff") << wolf.getCutoff() << " "+angstrom << endl
                  << pad(SUB,w,"Kappa") << wolf.getKappa() << " "+angstrom+"^-1" << endl;
                return o.str();
              }
          };
          
          class IonDipoleDamped : public IonDipole {
            private:
              string _brief() { return "Damped Ion-dipole"; }
              GaussianDamping gd;
            public:
              IonDipoleDamped(InputMap &in) : IonDipole(in),
              gd() {
                name+=" IonDipole";
              }

              template<class Tparticle>
                double operator()(const Tparticle &a, const Tparticle &b, const Point &r) const {
                  return _lB*gd.q2mu(a.charge*b.muscalar, b.mu, b.charge*a.muscal, a.mu, a.betaC, b.betaC, a.betaD, b.betaD,r);
                }

              template<class Tparticle>
                Point field(const Tparticle &p,const Tparticle &p0, const Point &r) const {
                  return Point(0,0,0);
                }

              string info(char w) {
                using namespace textio;
                std::ostringstream o;
                o << IonDipole::info(w) << endl;
                return o.str();
              }
          };

          class DipoleDipoleWolf : public DipoleDipole {
            private:
              string _brief() { return "Dipole-dipole Wolf"; }
              WolfBase wolf;
            public:
              DipoleDipoleWolf(InputMap &in) : DipoleDipole(in),
              wolf(in.get<double>("kappa", 1.8, "Kappa-damping"),
                  in.get<double>("dipdip_cutoff",in.get<double>("cuboid_len",pc::infty)/2)) {
                name+=" Wolf";
              }

              template<class Tparticle>
                double operator()(const Tparticle &a, const Tparticle &b, const Point &r) const {
                  return _lB*wolf.mu2mu(a.mu,b.mu, a.muscalar*b.muscalar, r);
                }

              /** @brief Dipole field at `r` due to dipole `p` 
               *  Gets returned in [e/Å] (\f$\beta eE \f$)
               */
              template<class Tparticle>
                Point field(const Tparticle &p,const Tparticle &p0, const Point &r) const {
                  return _lB*wolf.field(p,p0,r);
                }

              string info(char w) {
                using namespace textio;
                std::ostringstream o;
                o << DipoleDipole::info(w)
                  << pad(SUB,w,"Cutoff") << wolf.getCutoff() << " "+angstrom << endl
                  << pad(SUB,w,"Kappa") << wolf.getKappa() << " "+angstrom+"^-1" << endl;
                return o.str();
              }
          };

          class DipoleDipoleDamped : public DipoleDipole {
            private:
              string _brief() { return "Damped Dipole-dipole"; }
              GaussianDamping gd;
            public:
              DipoleDipoleDamped(InputMap &in) : DipoleDipole(in),
              gd() {
                name+=" DipoleDipole";
              }

              template<class Tparticle>
                double operator()(const Tparticle &a, const Tparticle &b, const Point &r) const {
                  return _lB*gd.mu2mu(a.charge*b.muscalar, b.mu, b.charge*a.muscal, a.mu, a.betaC, b.betaC, a.betaD, b.betaD,r);
                }

              template<class Tparticle>
                Point field(const Tparticle &p,const Tparticle &p0, const Point &r) const {
                  return _lB*gd.fieldMu2Mu(p,p0,r);
                }

              string info(char w) {
                using namespace textio;
                std::ostringstream o;
                o << DipoleDipole::info(w) << endl;
                return o.str();
              }
          };

          class DipoleDipoleWolfDamped : public DipoleDipole {
            private:
              string _brief() { return "Dipole-dipole Wolf Damped"; }
              WolfGaussianDamping wolf;
            public:
              DipoleDipoleWolfDamped(InputMap &in) : DipoleDipole(in),
              wolf(in.get<double>("kappa", 1.8, "Kappa-damping"),
                  in.get<double>("dipdip_cutoff",in.get<double>("cuboid_len",pc::infty)/2)) {
                name+=" Wolf Damped";
              }

              template<class Tparticle>
                double operator()(const Tparticle &a, const Tparticle &b, const Point &r) const {
                  return _lB*wolf.mu2mu(a.mu,b.mu, a.muscalar*b.muscalar, a.id-1, b.id-1, r);
                }

              /** @brief Dipole field at `r` due to dipole `p` 
               *  Gets returned in [e/Å] (\f$\beta eE \f$)
               */
              template<class Tparticle>
                Point field(const Tparticle &p,const Tparticle &p0, const Point &r) const {
                  return _lB*wolf.field(p,p.id,r);
                }

              string info(char w) {
                using namespace textio;
                std::ostringstream o;
                o << DipoleDipole::info(w)
                  << pad(SUB,w,"Cutoff") << wolf.getCutoff() << " "+angstrom << endl
                  << pad(SUB,w,"Kappa") << wolf.getKappa() << " "+angstrom+"^-1" << endl;
                return o.str();
              }
          };

          class IonQuadWolfDamped : public IonQuad {
            private:
              string _brief() { return "Ion-quad Wolf Damped"; }
              WolfGaussianDamping wolf;
            public:
              IonQuadWolfDamped(InputMap &in) : IonQuad(in),
              wolf(in.get<double>("kappa", 1.8, "Kappa-damping"),
                  in.get<double>("dipdip_cutoff",in.get<double>("cuboid_len",pc::infty)/2)) {
                name+=" Wolf Damped";
              }

              template<class Tparticle>
                double operator()(const Tparticle &a, const Tparticle &b, const Point &r) const {
                  return _lB*wolf.q2quad(a.charge, a.theta, b.charge, b.theta,a.id-1,b.id-1, r);
                }
                
              template<class Tparticle>
                Point field(const Tparticle &p,const Tparticle &p0, const Point &r) const {
                  return _lB*wolf.field(p,p0,r);
                }

              string info(char w) {
                using namespace textio;
                std::ostringstream o;
                o << IonQuad::info(w)
                  << pad(SUB,w,"Cutoff") << wolf.getCutoff() << " "+angstrom << endl
                  << pad(SUB,w,"Kappa") << wolf.getKappa() << " "+angstrom+"^-1" << endl;
                return o.str();
              }
          };


          class IonQuadWolf : public IonQuad {
            private:
              string _brief() { return "Ion-Quadrupole Wolf"; }
              WolfBase wolf;
            public:
              IonQuadWolf(InputMap &in) : IonQuad(in),
              wolf(in.get<double>("kappa", 1.8, "Kappa-damping"),
                  in.get<double>("dipdip_cutoff",in.get<double>("cuboid_len",pc::infty)/2)) {
                name+=" Wolf";
              }

              template<class Tparticle>
                double operator()(const Tparticle &a, const Tparticle &b, const Point &r) const {
                  return _lB*wolf.q2quad(a.charge, b.theta, b.charge, a.theta,r);
                }

              template<class Tparticle>
                Point field(const Tparticle &p,const Tparticle &p0, const Point &r) const {
                  return Point(0,0,0);
                }

              string info(char w) {
                using namespace textio;
                std::ostringstream o;
                o << IonQuad::info(w)
                  << pad(SUB,w,"Cutoff") << wolf.getCutoff() << " "+angstrom << endl
                  << pad(SUB,w,"Kappa") << wolf.getKappa() << " "+angstrom+"^-1" << endl;
                return o.str();
              }
          };
          
        }
        }
#endif

