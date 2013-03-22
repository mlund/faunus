#ifndef FAUNUS_MULTIPOLE
#define FAUNUS_MULTIPOLE

namespace Faunus {

  template<typename Tenergy>
    void getInducedDipoles(Tenergy &pot, Space &spc, const Eigen::MatrixXd &mu_per, double &limit) { 
      using namespace Eigen;
      int size = spc.p.size();
      double diLim = size*limit;   double diConv;
      MatrixXd Field(size,3);      Field.setZero();
      Point E(0,0,0);              Point mu_ind(0,0,0);    Point mu_err(0,0,0);

      do {  
        diConv = 0.0;
        Field = getField(pot,spc);
        for(int i = 0; i < size; i++) {
          E = Field.row(i);
          mu_ind = atom[spc.p[i].id].alphamatrix*E;
          mu_err = mu_ind - spc.p[i].mu*spc.p[i].muscalar + mu_per.row(i).transpose();
          diConv = diConv + mu_err.norm();
          mu_ind = mu_ind + mu_per.row(i).transpose();

          spc.p[i].muscalar = mu_ind.norm();
          spc.p[i].mu = mu_ind/spc.p[i].muscalar;
        }
      } while (diConv > diLim);
    }


  /**
   * @brief Returns dipole-dipole interaction which times the 
   * @brief bjerrum length in (Å) gives the energy in k_BT.
   *
   * @param mu1 Unit dipole moment vector 1
   * @param mu2 Unit dipole moment vector 2
   * @param mu1xmu2 Product of dipole scalars
   * @param r Vector R_{1->2}
   */
  template<class Tvec>
    double mu2mu(const Tvec &mu1, const Tvec &mu2, double mu1xmu2, const Tvec &r) {
      using namespace Eigen;
      double R1 = 1.0/r.norm();        //   R^-1 , Å^-1
      double R3 = R1*R1*R1;           //   R^-3 , Å^-3
      double R5 = 3.0*R3*R1*R1;       // 3*R^-5 , Å^-5
      Matrix3d T = R5*r*r.transpose() - R3*Matrix3d::Identity();  // T_{ab} = 3*r_a*r_b / |r|^5 , for a ne b , Å^-3
      // T_{aa} = 3*r_a^2 / |r|^5 - 1 / |r|^3    , Å^-3
      double W = -mu1.transpose()*T*mu2;                        // Summation of all T_{ab}*mu_a*mu_b       , Å^-3
      return W*mu1xmu2;                             // Scaling with |mu_a|*|mu_b|              , Å^-1 e^ 2
    }

  namespace Potential {
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

        Point field(const particle &p, const Point &r) const {
          double R1 = 0.0;          double R2 = 0.0;   double R3 = 0.0;                                                        
          Point E(0,0,0);                                                                         
          R1 = 1.0/r.norm();
          Point rn = r*R1;   R2 = R1*R1;   R3 = R2*R1;                                                       
          E = E + p.charge*R2*rn;                                             // From charges                        
          E = E + (3.0*p.mu.dot(rn)*rn - p.mu)*p.muscalar*R3;    // From dipoles                        
          return E;                                                                                                        

        }

        string info(char w) { return _brief(); }
    };

    class DipoleDipoleRF : public DipoleDipole {
      private:
        string _brief() { return "Reaction Field"; }
        double rc2,eps;
      public:
        DipoleDipoleRF(InputMap &in) : DipoleDipole(in) { 
          rc2 = in.get<double>("dipdip_cutoff",pc::infty);
          rc2 = rc2*rc2;
          eps = in.get<double>("epsilon_rf",80.);
          eps = _lB*(2*(eps-1)/(eps+1))/pow(rc2,1.5);
        }
        template<class Tparticle>
          double operator()(const Tparticle &a, const Tparticle &b, const Point &r) const {
            //cout << DipoleDipole::operator()(a,b,r) << " " << eps*a.mu.dot(b.mu)*a.muscalar*b.muscalar << "\n";
            if (r.squaredNorm() < rc2)
              return DipoleDipole::operator()(a,b,r) - eps*a.mu.dot(b.mu)*a.muscalar*b.muscalar;
            return 0.0;
          }

        string info(char w) { return _brief(); }
    };
  }
}
#endif

