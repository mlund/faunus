#ifndef FAUNUS_MULTIPOLE
#define FAUNUS_MULTIPOLE

namespace Faunus {

  template<typename Tenergy>
    void getInducedDipoles(Tenergy &pot, Space &spc, const Eigen::MatrixXd &mu_per, const double &limit) { 
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
          mu_ind = spc.p[i].alpha*E;
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
   * @param a DipoleParticle 1
   * @param b DipoleParticle 2
   * @param r Vector R_{1->2}
   */
  double mu2mu(const particle &a, const particle &b, const Point &r) {
    using namespace Eigen;
    double R1 = 1.0/r.norm();        //   R^-1 , Å^-1
    double R3 = -R1*R1*R1;           //  -R^-3 , Å^-3
    double R5 = -3.0*R3*R1*R1;       // 3*R^-5 , Å^-5

    Matrix3d T = R5*r*r.transpose() + R3*Matrix3d::Identity();  // T_{ab} = 3*r_a*r_b / |r|^5 , for a ne b , Å^-3
    // T_{aa} = 3*r_a^2 / |r|^5 - 1 / |r|^3    , Å^-3
    double W = -a.mu.transpose()*T*b.mu;                        // Summation of all T_{ab}*mu_a*mu_b       , Å^-3
    return W*a.muscalar*b.muscalar;                             // Scaling with |mu_a|*|mu_b|              , Å^-1 e^ 2
  }

  namespace Potential {
    class DipDip : public PairPotentialBase {
      private:
        string _brief() { return "Dipole-dipole"; }
        double _lB;
      public:
        DipDip(InputMap &in) {
          _lB=in.get<double>("bjerrumlength", 7.1);
        }
        double operator()(const particle &a, const particle &b, const Point &r) const {
          return _lB*mu2mu(a,b,r);
        }
        double operator()(const particle &a, const particle &b, double r2) const {
          assert(!"Should not reach here!");
          return pc::infty;
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
  }
}
#endif

