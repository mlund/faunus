#ifndef FAUNUS_MULTIPOLE_H
#define FAUNUS_MULTIPOLE_H

namespace Faunus {
  
  template<class Tvec>
    double nemoRep(double aab, double bab, double cab, double dab, double eab, double fab, const Tvec &r) {
      double r2i = 1/r.squaredNorm();
      double r6i = r2i*r2i*r2i;
      double uexp1  = eab*pow(2.71828,-std::min(80.0,fab/sqrt(r2i)));
      double uexp2  = dab*r6i*r6i*r2i;
      double udis1  =-cab*r6i;
      double udis2  = aab*pow(2.71828,-std::min(80.0,bab/sqrt(r2i)));
      double ur = (uexp1 + uexp2  + udis1 + udis2);
      // (q1*q2*lB/r)*kT*N_AV/1000 = kJ/mol
      //return (1000/(pc::Nav*pc::kT()))*ur;
      return ur;
    }
  
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
   * @param muA Unit dipole moment vector of particle A
   * @param muB Unit dipole moment vector of particle B
   * @param muAxmuB Product of dipole scalars
   * @param r Vector \f$ r_{AB} \f$
   */
  template<class Tvec>
    double mu2mu(const Tvec &muA, const Tvec &muB, double muAxmuB, const Tvec &r) {
      double R2 = 1/r.squaredNorm();
      double R1 = sqrt(R2);
      double R3 = R1*R2;
      double R5 = R3*R2;
      //Eigen::Matrix3d T = 3*R5*r*r.transpose() - R3*Matrix3d::Identity();
      //double W = -muA.transpose()*T*muB;                       // Buckingham    Å^-3
      double W = muA.dot(muB)*R3-3*muA.dot(r)*muB.dot(r)*R5; // J&K
      return W*muAxmuB;  // e^2 Å^2 Å ^-3 = e^2 /A
    }

  /**
   * @brief Returns ion-quadrupole interaction
   */
  template<class Tvec, class Tmat>
    double q2quad(double q, const Tmat &quad, const Tvec &r) {
      double R2 = 1/r.squaredNorm();
      double R1 = sqrt(R2);
      double R3 = R1*R2;
      double R5 = R3*R2;
      double W = r.transpose()*quad*r;
      W = W*R5  - quad.trace()*(R3/3); // e / Å
      return q*W; // e^2 / Å
    }

  namespace Potential {
    
    /**
     * @brief Ion-dipole interaction, 
     *
     * More info...
     */
    class NemoRepulsion : public PairPotentialBase {
      private:
        string _brief() { return "NemoRepulsion"; }
      protected:
        typedef Eigen::VectorXd Tvec;
        typedef opair<int> Tpair;
        std::map<Tpair,Tvec> map;

        double _lB, epsilon_r, rc2, eps;
        Eigen::MatrixXd pab;
        particle::Tid HW,OW; // particle ID

      public:
        NemoRepulsion(InputMap &in) {
          name="Nemo repulsion";
          pc::setT ( in.get<double>("temperature", 298.15, "Absolute temperature (K)") );
          epsilon_r = in.get<double>("epsilon_r",80., "Dielectric constant");
          _lB=pc::lB( epsilon_r );
          rc2 = pow(in.get<double>("dipdip_cutoff",pc::infty), 2);
          eps = in.get<double>("epsilon_rf",80.);
          eps = _lB*(2*(eps-1)/(eps+1))/pow(rc2,1.5);
          pab.resize(3,6);
          pab << 28063.645684, 3.816817, 816.621758, 13824.695518, 518387.587215, 4.194060,
              1872.664555,  3.997893, 214.610684, 8857.398499,  11055.163586,  3.534167,
              597.625692,   4.197006, 16.5577400, 0.000000,     800.775380,    2.677176;
          pab = pab*(1000/(pc::Nav*pc::kT()));
          pab.col(1) = pab.col(1)/(1000/(pc::Nav*pc::kT()));
          pab.col(5) = pab.col(5)/(1000/(pc::Nav*pc::kT()));
          HW = atom["HW"].id;
          OW = atom["OW"].id;
        }

        /**
         * @brief NemoRepulsion
         * @param a Dipole particle A
         * @param b Dipole particle B
         * @param r Direction \f$ r_A - r_B \f$  
         */
        template<class Tparticle> // q2mu(1->2,r) + q2mu(2->1,-r) = q2mu(1->2,r) - q2mu(2->1,r)
          double operator()(const Tparticle &a, const Tparticle &b, const Point &r) const {
            int temp = 1;
            if(a.id == HW && b.id == HW) {
              temp = 2;
            } else if(a.id == OW && b.id == OW) {
              temp = 0;
            }
            //double iondipole_e = _lB*(q2mu(a.charge*b.muscalar,b.mu,r) - q2mu(b.charge*a.muscalar,a.mu,r));
            //double dipdip_e = _lB*mu2mu(a.mu, b.mu, a.muscalar*b.muscalar, r);
            double nemoRep_e = nemoRep(pab(temp,0),pab(temp,1),pab(temp,2),pab(temp,3),pab(temp,4),pab(temp,5), r);
            //double ionquad_e = _lB*(q2quad(a.charge, b.theta,r)+q2quad(b.charge, a.theta,r));
            //if (r.squaredNorm() < rc2)
            //  dipdip_e = dipdip_e - eps*a.mu.dot(b.mu)*a.muscalar*b.muscalar;
            return (nemoRep_e);
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
        string _brief() {
          std::ostringstream o;
          o << "Dipole-dipole, lB=" << _lB << textio::_angstrom;
          return o.str();          
        }
      protected:
        double _lB;
        double convert;
      public:
        DipoleDipole(InputMap &in) {
          pc::setT ( in.get<double>("temperature", 298.15,
                "Absolute temperature (K)") );
          double epsilon_r = in.get<double>("epsilon_r",80.,
              "Dielectric constant");
          _lB = pc::lB(epsilon_r);
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
              return DipoleDipole::operator()(a,b,r)
                - eps*a.mu.dot(b.mu)*a.muscalar*b.muscalar;
            return 0;
          }

        void updateDiel(double er) {
          eps = _lB*(2*(er-1)/(er+1))/pow(rc2,1.5);
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
            return _lB*(q2quad(a.charge, b.theta,r)+q2quad(b.charge, a.theta,r));
          }

        string info(char w) { return _brief(); }
    };
  }
}
#endif

