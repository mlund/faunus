#include <faunus/faunus.h>

using namespace Faunus;
using namespace std;
using namespace Eigen;

// See include/faunus/potential.h

namespace Faunus {
  namespace Nemo {

    double getSubTensor(Matrix3d M, const Point &r, const VectorXi &v) {
      if (v[1] == 1) {
        M = M/r.x();
        return M(v[2],v[3]);
      } else if (v[1] == 2) {
        M = M/r.y();
        return M(v[2],v[3]);
      } else if (v[1] == 3) {
        M = M/r.z();
        return M(v[2],v[3]);
      } else {
        return 0.0;
      }	
    }

    double getTensor(const Point &r, const VectorXi &v) {
      double R = r.norm();
      double R2 = 1.0/(R*R);
      double T = 0.0;
      Vector3d V(R/r.x(), R/r.y(), R/r.z());
      Matrix3d M; M << R/(r.x()*r.x()), R/(r.x()*r.y()), R/(r.x()*r.z()), R/(r.y()*r.x()),
               R/(r.y()*r.y()), R/(r.y()*r.z()), R/(r.z()*r.x()), R/(r.z()*r.y()),
               R/(r.z()*r.z());
      int size = v.size();
      switch (size) {
        case 1:
          T = -R2*V[v[0]];
          break;
        case 2:
          T = (2/R)*R2*V[v[0]]*V[v[1]] - R2*M(v[0],v[1]);	
          break;
        case 3:
          T = (-6.0*R2*R2)*V[v[0]]*V[v[1]]*V[v[2]];
          T = T + (2.0/R)*R2*M(int(v[1]),int(v[2]))*V[v[0]];
          T = T + (2.0/R)*R2*M(int(v[2]),int(v[0]))*V[v[1]];
          T = T + (2.0/R)*R2*M(int(v[0]),int(v[1]))*V[v[2]];
          T = T + (-R2*getSubTensor(M,r,v));
          break;
      }
      return T;
    }

    double q2q(const particle &a, const particle &b, const Point &r) {
      return a.charge*b.charge/r.norm();
    }

    double q2mu(const particle &a, const particle &b, const Point &r) {
      double W = 0.0;
      VectorXi V(0);
      for(int i = 0; i < 3; i++) {
        V(0) = i;
        W = W + getTensor(r,V)*(a.charge*b.mu(i)*b.muscalar + b.charge*a.mu(i)*a.muscalar);
      }
      return W;
    }

    double mu2mu(const particle &a, const particle &b, const Point &r) {
      double W = 0.0;
      Vector2i V(0,0);
      double sca = double(1/3);
      for(int i = 0; i < 3; i++) {
        V(0) = i;
        for(int j = i; j < 3; j++) {
          V(1) = j;
          W = W + getTensor(r,V)*(sca*a.charge*b.theta(i,j) - a.mu(i)*a.muscalar*b.mu(j)*b.muscalar + sca*a.theta(i,j)*b.charge);
        }
      }	
      return W;
    }

    double mu2theta(const particle &a, const particle &b, const Point &r) {
      double W = 0.0;
      Vector3i V(0,0,0);
      double sca = double(1/3);
      for(int i = 0; i < 3; i++) {
        V(0) = i;
        for(int j = i; j < 3; j++) {
          V(1) = j;
          for(int k = j; k < 3; k++) {
            V(2) = k;
            W = W + sca*getTensor(r,V)*(a.theta(i,j)*b.mu(k)*b.muscalar - a.mu(i)*a.muscalar*b.theta(j,k));
          }
        }
      } 
      return W;
    }
  }//namespace nemo
}//namespace faunus

namespace Faunus {
  namespace Potential {

    class DipoleDipole : public PairPotentialBase {
      private:
        string _brief() { return "Dipole-dipole"; }
        double _lB;
      public:
        DipoleDipole(InputMap &in) {
          _lB=in.get<double>("bjerrumlength", 7.1);
        }
        double operator()(const particle &a, const particle &b, const Point &r) const {
          return _lB*Faunus::Nemo::mu2mu(a,b,r);
        }
        string info(char w) { return _brief(); }
    };

    typedef CombinedPairPotential<HardSphere,DipoleDipole> THSDipDip;

  }//namespace
}//namespace

int main() {
  // particle a;
  // a.mu.x()=1.0;
  // Point b(3,4,1);
  // VectorXi v; v << 1,1,1;
  // double T = Nemo::getTensor(b,v);
  // cout << T << endl;
  typedef Potential::THSDipDip Tpair; // d(mcp);

  typedef Geometry::Cuboid Tgeo;                        // select simulation geometry and pair potential
  // typedef Potential::CoulombLJ Tpair;
  ::atom.includefile("minimal.json");                 // load atom properties
  InputMap in("minimal.input");                       // open parameter file for user input
  Energy::NonbondedVector<Tpair,Tgeo> pot(in);              // create Hamiltonian, non-bonded only
  Space spc( pot.getGeometry() );                     // create simulation space, particles etc.
  GroupAtomic salt(spc, in);                          // group for salt particles
  Move::AtomicTranslation mv(in, pot, spc);           // particle move class
  mv.setGroup(salt);                                  // tells move class to act on salt group
  mv.move(1e3);                                       // move salt randomly 100000 times
  std::cout << spc.info() << pot.info() << mv.info(); // final information
  std::cout << spc.p[0].muscalar << endl;
}


