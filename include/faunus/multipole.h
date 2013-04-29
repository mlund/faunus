#ifndef FAUNUS_MULTIPOLE_H
#define FAUNUS_MULTIPOLE_H

namespace Faunus {

  /**
   * @brief Returns ion-dipole interaction, Needs to be checked!
   * @param QxMu Product of ion charge and dipole scalar
   * @param mu Unit dipole moment vector
   * @param r Direction \f$ r_Mu - r_Q \f$  
   *
   */
  template<class Tvec>
    double q2mu(double QxMu, const Tvec &mu, const Tvec &r) {
      double R2 = 1/r.squaredNorm();
      double R1 = sqrt(R2);
      double R3 = R1*R2;
      double W = QxMu*mu.dot(r)*R3;     
      return W;  // Beware of r_Mu - r_Q = -r according to Israelachvili p.36, i.e. minus becomes plus
    }
  
  /**
   * @brief Returns dipole-dipole interaction
   *
   * @param mu1 Unit dipole moment vector of particle A
   * @param mu2 Unit dipole moment vector of particle B
   * @param mu1xmu2 Product of dipole scalars
   * @param r Vector \f$ r_{AB} \f$
   *
   */
  template<class Tvec>
    double mu2mu(const Tvec &muA, const Tvec &muB, double muAxmuB, const Tvec &r) {
      double R2 = 1/r.squaredNorm();
      double R1 = sqrt(R2);
      double R3 = R1*R2;
      double R5 = R3*R2;
      Eigen::Matrix3d T = 3*R5*r*r.transpose() - R3*Matrix3d::Identity();
      double W = -muA.transpose()*T*muB;                       // Buckingham
      //double W = mu1.dot(mu2)*R3-3*mu1.dot(r)*mu2.dot(r)*R5; // J&K
      return W*muAxmuB;
    }
    
  /**
   * @brief Returns ion-quadrupole interaction
   *
   * @param mu1 Charge of particle 1
   * @param mu2 Quadrupole matrix of particle 2
   * @param r Vector \f$ R_{ij} \f$
   *
   */
  template<class Tvec, class Tmat>
    double q2quad(double q, const Tmat &quad, const Tvec &r) {
      double R2 = 1/r.squaredNorm();
      double R1 = sqrt(R2);
      double R3 = R1*R2;
      double R5 = R3*R2;
      double W = r.transpose()*quad*r;
      W = W*R5  - quad.trace()*(R3/3);
      return q*W;
    }
    
  namespace Potential {

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
            return _lB*(q2mu(a.charge*b.muscalar,b.mu,r) - q2mu(b.charge*a.muscalar,a.mu,r)); 
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
        string _brief() { return "Dipole-dipole"; }
      protected:
        double _lB;
	double convert;
    public:
        DipoleDipole(InputMap &in) {
          pc::setT ( in.get<double>("temperature", 298.15, "Absolute temperature (K)") );
          double epsilon_r = in.get<double>("epsilon_r",80., "Dielectric constant");
          _lB=pc::lB( epsilon_r );
	  convert = _lB*pc::kT()/(pc::e*pc::e);
        }
        template<class Tparticle>
          double operator()(const Tparticle &a, const Tparticle &b, const Point &r) const {
            return _lB*mu2mu(a.mu, b.mu, a.muscalar*b.muscalar, r);
          }

        /** @brief Dipole field at `r` due to dipole `p` 
         *  Gets returned in [e/Å] (\f$\beta eE \f$)
         */
        template<class Tparticle>
          Point field (const Tparticle &p, const Point &r) const {
            double R2 = 1.0/r.squaredNorm();
            double R1 = sqrt(R2);
            Point r_n = r*R1;
            return ((3.0*p.mu.dot(r_n)*r_n - p.mu)*R2*R1)*p.muscalar*_lB; // \beta e E
          }
        string info(char w) { return _brief(); }
    };

    /**
     * @brief Dipole-dipole interaction w. spherical cutoff and reaction field
     *
     * More info...
     */
    class DipoleDipoleRF : public DipoleDipole {
      private:
        string _brief() { return "Dipole-dipole (RF)"; }
        double rc2,eps;
      public:
        DipoleDipoleRF(InputMap &in) : DipoleDipole(in) { 
          rc2 = pow(in.get<double>("dipdip_cutoff",pc::infty), 2);
          eps = in.get<double>("epsilon_rf",80.);
          eps = _lB*(2*(eps-1)/(eps+1))/pow(rc2,1.5);
        }
        template<class Tparticle>
          double operator()(const Tparticle &a, const Tparticle &b, const Point &r) const {
            if (r.squaredNorm() < rc2)
              return DipoleDipole::operator()(a,b,r) - eps*a.mu.dot(b.mu)*a.muscalar*b.muscalar;
            return 0.0;
          }

        string info(char w) { return _brief(); }
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
          pc::setT ( in.get<double>("temperature", 298.15, "Absolute temperature (K)") );
          double epsilon_r = in.get<double>("epsilon_r",80., "Dielectric constant");
          _lB=pc::lB( epsilon_r );
        }
        template<class Tparticle>
          double operator()(const Tparticle &a, const Tparticle &b, const Point &r) const {
            return _lB*q2quad(a.charge, b.theta,r);
          }

        string info(char w) { return _brief(); }
    };
    
  }
}
#endif

