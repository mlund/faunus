#include <faunus/faunus.h>
#include <faunus/physconst.h>
#include <math.h>
#include <faunus/potentials.h>
#include <faunus/geometry.h>

using namespace Faunus;
using namespace std;
using namespace Eigen;

// See include/faunus/potential.h

namespace Faunus {
  namespace Nemo {

    /**
     * @brief Returns charge-dipole interaction which times the
     * @brief bjerrum length in (Å) gives energy in k_BT.
     *
     * @param a DipoleParticle 1
     * @param b DipoleParticle 2
     * @param r Vector R_{1->2}
     */
    double q2mu(const particle &a, const particle &b, const Point &r) {
      double R1 = 1.0/r.norm();        //   R^-1 , Å^-1
      double R3 = -R1*R1*R1;           //  -R^-3 , Å^-3

      double W1 = r.transpose()*b.mu;  //          Å
      double W2 = r.transpose()*a.mu;  //          Å
      double W = (W1*b.muscalar*a.charge - W2*a.muscalar*b.charge)*R3; 
      return W;                        //          Å^-1 e^ 2
    }

    /**
     * @brief Returns dipole-dipole interaction which times the 
     * @brief bjerrum length in (Å) gives the energy in k_BT.
     *
     * @param a DipoleParticle 1
     * @param b DipoleParticle 2
     * @param r Vector R_{1->2}
     */

    /*
    double mu2mu(const particle &a, const particle &b, const Point &r) {
      double R1 = 1.0/r.norm();        //   R^-1 , Å^-1
      double R3 = -R1*R1*R1;           //  -R^-3 , Å^-3
      double R5 = -3.0*R3*R1*R1;       // 3*R^-5 , Å^-5

      Matrix3d T = R5*r*r.transpose() + R3*Matrix3d::Identity();  // T_{ab} = 3*r_a*r_b / |r|^5 , for a ne b , Å^-3
                                                                  // T_{aa} = 3*r_a^2 / |r|^5 - 1 / |r|^3    , Å^-3
      double W = -a.mu.transpose()*T*b.mu;                        // Summation of all T_{ab}*mu_a*mu_b       , Å^-3
      return W*a.muscalar*b.muscalar;                             // Scaling with |mu_a|*|mu_b|              , Å^-1 e^ 2
    } 
    */

    /*
    double q2theta(const particle &a, const particle &b, const Point &r) {
      double R1 = 1.0/r.norm();        //   R^-1 , Å^-1
      double R3 = -R1*R1*R1;           //  -R^-3 , Å^-3
      double R5 = -3.0*R3*R1*R1;       // 3*R^-5 , Å^-5

      Matrix3d T = R5*r*r.transpose() + R3*Matrix3d::Identity();  // T_{ab} = 3*r_a*r_b / |r|^5 , for a ne b , Å^-3
      double W1 = (a.theta.transpose()*T).trace()*b.charge/3.0;   // e^2 Å^-1   
      double W2 = (b.theta.transpose()*T).trace()*a.charge/3.0;      
      double W = W1 + W2;
      return W;  //
    }
    */
  }//namespace nemo
}//namespace faunus
/*
namespace Faunus {
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
          return _lB*Faunus::Nemo::mu2mu(a,b,r);
        }
        string info(char w) { return _brief(); }
    };

    //typedef CombinedPairPotential<HardSphere,DipoleDipole> THSDipDip;

  }//namespace
}//namespace
*/
double getKeesom(Space &spc, double R, double lB) {
  double keesom = -pow(spc.p[0].muscalar,2)*pow(spc.p[1].muscalar,2)*lB*lB/3;
  return keesom/pow(R,6);
}

int main() {
  //typedef Potential::CombinedPairPotential<Potential::HardSphere, Potential::DipoleDipole> Tpair;
  typedef Potential::DipoleDipoleHS Tpair;                 // Select potential
  //typedef Potential::HardSphere Tpair;                 // Select potential
  typedef Geometry::Cuboid Tgeo;                      // Select simulation geometry and pair potential
  ::atom.includefile("minimal.json");                 // Load atom properties
  InputMap in("minimal.input");                       // Open parameter file for user input
  Energy::NonbondedVector<Tpair,Tgeo> pot(in);        // Create Hamiltonian, non-bonded only
  Space spc( pot.getGeometry() );                     // Create simulation space, particles etc.
  GroupAtomic salt(spc, in);                          // Group for salt particles

  // Initialization of particles
  spc.p[0].z() = -10;                    // In Å
  spc.p[1].z() = 10;                     // In Å
  for(size_t i = 0; i < spc.p.size(); i++) {
    spc.p[i].x() = 0;
    spc.p[i].y() = 0;
    //spc.p[i].theta << -4.9556, -0.0000, -0.0000,-0.0000,-7.7320,0.0000, -0.0000, 0.0000, -6.0094; // In Debye*Å
    spc.trial[i] = spc.p[i];
  }

  Analysis::Table2D<float,unsigned long int> rdf(0.2);
  Move::AtomicTranslation mv(in, pot, spc);           // particle move class
  Move::AtomicRotation rot(in, pot, spc);           // particle move class
  mv.setGroup(salt);                                  // tells move class to act on salt group
  rot.setGroup(salt);                                  // tells move class to act on salt group
  mv.dir=Point(0,0,1);

  int N = 1000000; // Number of random displacments of dipole moment per fixed distance

  Potential::DipoleDipole dd(in);

  for(int i = 0; i < N; i++) {
    if (slp_global()>0.5)
      mv.move(salt.size());            // move salt randomly 100000 times
    else
      rot.move(salt.size());            // move salt randomly 100000 times

    double R=spc.geo->dist(spc.p[0],spc.p[1]);
    rdf(R)++;
  }

  rdf.save("rdf_nemo.dat");
  std::cout << spc.info() << pot.info() << mv.info() << rot.info(); // final information
  return 1;
}


MatrixXd getField(Energy::NonbondedVector<Potential::DipoleDipoleHS,Geometry::Cuboid> &pot, Space &spc) { 
  int size = spc.p.size();
  double R1 = 0.0;          double R2 = 0.0;   double R3 = 0.0;
  Point E(0,0,0);           Point r(0,0,0);
  MatrixXd field(size,3);   field.setZero();

  for(int I =0; I < size; I ++) {
    for(int i = 0; i < size; i ++) {
      if(i == I) continue;
      r = pot.getGeometry().vdist(spc.p[i],spc.p[I]);
      R1 = 1.0/r.norm();   r = r*R1;   R2 = R1*R1;   R3 = R2*R1;
      E = E + spc.p[i].charge*R2*r;                                             // From charges
      E = E + (3.0*spc.p[i].mu.dot(r)*r - spc.p[i].mu)*spc.p[i].muscalar*R3;    // From dipoles
      field.row(I) = E;
    }
    E.setZero();
  }
  return field;
}
/*
void getInducedDipoles(Energy::NonbondedVector<Potential::DipoleDipoleHS,Geometry::Cuboid> &pot, Space &spc, const MatrixXd &mu_per, const double &limit) { 
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
*/

int main1() {
  //typedef Potential::CombinedPairPotential<Potential::HardSphere, Potential::DipoleDipole> Tpair;
  typedef Potential::DipoleDipoleHS Tpair;                 // Select potential
  typedef Geometry::Cuboid Tgeo;                      // Select simulation geometry and pair potential
  ::atom.includefile("minimal.json");                 // Load atom properties
  InputMap in("minimal.input");                       // Open parameter file for user input
  Energy::NonbondedVector<Tpair,Tgeo> pot(in);        // Create Hamiltonian, non-bonded only
  Space spc( pot.getGeometry() );                     // Create simulation space, particles etc.
  GroupAtomic salt(spc, in);                          // Group for salt particles

  for(size_t i = 0; i < spc.p.size(); i++) {
    spc.p[i].muscalar = 1.85*0.20819434; // In eÅ for water
    //spc.p[i].theta << -4.9556, -0.0000, -0.0000,-0.0000,-7.7320,0.0000, -0.0000, 0.0000, -6.0094; // In Debye*Å
    //spc.p[i].alpha << 7.14264826, 0.00014943, -0.00003816, 0.00014943,5.64937691,0.00035990, -0.00003816, 0.00035990, 6.03919669; 
    //spc.p[i].alpha << 2.6,0,0,0,2.6,0,0,0,2.6; 
    //spc.p[i].alpha = spc.p[i].alpha;
    spc.p[i].x() = 0;
    spc.p[i].y() = 0;
    spc.trial[i] = spc.p[i];
  }
  spc.p[0].z() = 0;
  spc.p[1].z() = 10;
  //spc.p[0].mu = Point(0,0,0);
  //spc.p[1].mu = Point(0,0,0);
  //spc.p[0].muscalar = 0;
  //spc.p[1].muscalar = 0;

  //double T = 298.0;
  double eps_r = 78.7;
  //double beta = 1.0/(pc::kB*T);                               // J^-1
  //double lB = beta*pow(pc::e,2)/(4*pc::pi*pc::e0*eps_r); // meter

  double limit = 1e-8;
  int size = spc.p.size();
  MatrixXd mu_per(size,3); mu_per.setZero();
  for(int i = 0; i < size; i++) {
    mu_per.row(i) = spc.p[i].mu*spc.p[i].muscalar;
  }
 
  //getInducedDipoles(pot, spc, mu_per, limit);
}
