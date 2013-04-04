#ifndef FAUNUS_MULTIPOLE_H
#define FAUNUS_MULTIPOLE_H

namespace Faunus {

  /**
   * @brief Replaces dipole moment with permanent dipole moment plus induced dipole moment
   *
   * @param pot The potential including geometry
   * @param spc The space including the particles
   * @param mu_per Matrix with permanent dipole moments of all particles
   * @param limit Limit of sum of dipole-changes when iteration should stop
   */
  template<typename Tenergy>
    void getInducedDipoles(Tenergy &pot, Space &spc, const Eigen::MatrixXd &mu_per, double &limit) { 
      int size = spc.p.size();
      double diLim = size*limit;
      double diConv;
      Eigen::MatrixXd Field(size,3);
      Field.setZero();
      Point E(0,0,0);
      Point mu_ind(0,0,0);
      Point mu_err(0,0,0);

      do {  
        diConv = 0.0;
        Field = getField(pot,spc);
        for(int i = 0; i < size; i++) {
          E = Field.col(i) + ExternalField;
          //mu_ind = atom[spc.p[i].id].alphamatrix*E;
          mu_ind = 0.0001*E;
          mu_err = mu_ind - spc.p[i].mu*spc.p[i].muscalar + mu_per.col(i);
          diConv = diConv + mu_err.norm();
          mu_ind = mu_ind + mu_per.col(i);

          spc.p[i].muscalar = mu_ind.norm();
          spc.p[i].mu = mu_ind/spc.p[i].muscalar;
        }
      } while (diConv > diLim);
    }

  /**
   * @brief Returns the dielectric constant outside the cutoff limit. Only hold when using PBC and $\epsilon_{sur} = \epsilon$,
   * @brief [Neumann, M. (1983) Mol. Phys., 50, 841-858].
   *
   * @param pot The potential including geometry
   * @param spc The space including the particles
   * @param cutoff The cutoff of the reaction field
   */
  template<typename Tenergy>
    double getDielectricConstant(Tenergy &pot, Space &spc, double cutoff) { 
      using namespace Eigen;
      double M2 = 0;
      Point origin(0,0,0);
      double N = 0;

      for(unsigned int i = 0; i < spc.p.size(); i ++) {
        if( pot.getGeometry().dist(spc.p[i],origin) < cutoff) {
          M2 = M2 + spc.p[i].muscalar*spc.p[i].muscalar;
          N = N + 1;
        }
      }
      M2 = M2*pow(10,-20)/N;
      double Q = 0.25 + M2*pc::pi*spc.p.size()/(pot.getGeometry().getVolume()*pc::kB*pc::T());
      return ( Q + std::sqrt(Q*Q+0.5) );
    }

  /**
   * @brief Returns the electric field on all particles in a matrix
   *
   * @param pot The potential including geometry
   * @param spc The space including the particles
   */
  template<typename Tenergy>
    MatrixXd getField(Tenergy &pot, const Space &spc) {
      int size = spc.p.size();
      double R1 = 0.0;          double R2 = 0.0;   double R3 = 0.0;
      Point E(0,0,0);           Point r(0,0,0);
      MatrixXd field(3,size);   field.setZero();

      for(int I =0; I < size; I ++) {
        for(int i = 0; i < I; i ++) {
          r = pot.getGeometry().vdist(spc.p[i],spc.p[I]);
          R1 = 1.0/r.norm();   r = r*R1;   R2 = R1*R1;   R3 = R2*R1;
          E = E + spc.p[i].charge*R2*r;                                             // From charges
          E = E + (3.0*spc.p[i].mu.dot(r)*r - spc.p[i].mu)*spc.p[i].muscalar*R3;    // From dipoles
          field.col(I) = E;
        }
        for(int i = I + 1; i < size; i++) {
          r = pot.getGeometry().vdist(spc.p[i],spc.p[I]);
          R1 = 1.0/r.norm();   r = r*R1;   R2 = R1*R1;   R3 = R2*R1;
          E = E + spc.p[i].charge*R2*r;                                             // From charges
          E = E + (3.0*spc.p[i].mu.dot(r)*r - spc.p[i].mu)*spc.p[i].muscalar*R3;    // From dipoles
          field.col(I) = E;
        }
        E.setZero();
      }
      return field;
    }

  /**
   * @brief Returns dipole-dipole interaction which times the 
   *
   * @param mu1 Unit dipole moment vector 1
   * @param mu2 Unit dipole moment vector 2
   * @param mu1xmu2 Product of dipole scalars
   * @param r Vector R_{1->2}
   *
   * T_{ab} = 3*r_a*r_b / |r|^5 , for a ne b
   * T_{aa} = 3*r_a^2 / |r|^5 - 1 / |r|^3
   * Summation of all T_{ab}*mu_a*mu_b
   *
   */
  template<class Tvec>
    double mu2mu(const Tvec &mu1, const Tvec &mu2, double mu1xmu2, const Tvec &r) {
      double R2 = 1/r.squaredNorm();
      double R1 = sqrt(R2);
      double R3 = R1*R2;
      double R5 = R3*R2;
      Eigen::Matrix3d T = 3*R5*r*r.transpose() - R3*Matrix3d::Identity();
      double W = -mu1.transpose()*T*mu2;
      return W*mu1xmu2;
    }
/*
  namespace Induced {
    class InduceDipoles {
      private:
        string _brief() { return "Induced dipoles"; }
      public:
        InducedDipoles(InputMap &in, const Space &spc) {
          double epsilon_r = in.get<double>("epsilon_r",80., "Dielectric constant");
          Eigen::MatrixXd mu_per(3,spc.p.size());
          for(unsigned int i = 0; i < spc.p.size(); i++) {
            mu_per.row(i) = spc.p[i].mu*spc.p[i].muscalar;
          }
        }
        template<class Tenergy>
          void operator()(const Tenergy &pot, const Space &spc, double limit, Point ExternalField) const { 
            induceDipoles(pot,spc,mu_per,limit,ExternalField);
            return;
          }


        string info(char w) { return _brief(); }
    };
  }
*/
  namespace Potential {

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
      public:
        DipoleDipole(InputMap &in) {
          pc::setT ( in.get<double>("temperature", 298.15, "Absolute temperature (K)") );
          double epsilon_r = in.get<double>("epsilon_r",80., "Dielectric constant");
          _lB=pc::lB( epsilon_r );
        }
        template<class Tparticle>
          double operator()(const Tparticle &a, const Tparticle &b, const Point &r) const {
            return _lB*mu2mu(a.mu,b.mu,a.muscalar*b.muscalar,r);
          }

        /** @brief Dipole field at `r` due to dipole `p` */
        template<class Tparticle>
        Point field (const Tparticle &p, const Point &r) const {
          double R2 = 1.0/r.squaredNorm();
          double R1 = sqrt(R2);
          Point r_n = r*R1;
          return _lB*((3.0*p.mu.dot(r_n)*r_n - p.mu)*p.muscalar*R2*R1);
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
  }
}
#endif

