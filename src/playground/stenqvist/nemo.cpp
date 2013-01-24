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
          T = -R2*V[v[1]];
          break;
        case 2:
          T = (2/R)*R2*V[v[1]]*V[v[2]] - R2*M(v[1],v[2]);	
          break;
        case 3:
          T = (-6.0*R2*R2)*V[v[1]]*V[v[2]]*V[v[3]];
          T = T + (2.0/R)*R2*M(int(v[2]),int(v[3]))*V[v[1]];
          T = T + (2.0/R)*R2*M(int(v[3]),int(v[1]))*V[v[1]];
          T = T + (2.0/R)*R2*M(int(v[1]),int(v[2]))*V[v[1]];
          T = T + (-R2*getSubTensor(M,r,v));
          break;
      }
      return T;
    }

    double IonDipole(const particle &a, const particle &b, const Point &r) {

      return 0.0;
    }

    double mu2mu(const particle &a, const particle &b, const Point &r) {

      //Point r = a.coord-b.coord;
      double W = 0.0;
      for(int i = 0; i < 3; i++) {

      }	
      return 0.0;
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
          return 0;
        }
        string info(char w) { return _brief(); }
    };

    typedef CombinedPairPotential<HardSphere,DipoleDipole> THSDipDip;

  }//namespace
}//namespace

int main() {
  InputMap mcp("nemo.input");
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
  mv.move(1e5);                                       // move salt randomly 100000 times
  std::cout << spc.info() << pot.info() << mv.info(); // final information
}


