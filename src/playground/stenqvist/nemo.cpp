#include <faunus/faunus.h>

// Manual: http://faunus.sourceforge.net/doxyhtml/index.html

using namespace Faunus;
using namespace std;

class DipoleParticle : public PointParticle {
  public:
    double mu_s;            //!< Dipole moment scalar
    Eigen::Vector3d mu;     //!< Dipole moment unit vector
    Eigen::Matrix3d alpha;  //!< Polarization tensor

    inline DipoleParticle() {
      mu_s=0;
    }

    template<typename OtherDerived>                                   
      DipoleParticle(const Eigen::MatrixBase<OtherDerived>& other) : Tvec(other) {}  

    template<typename OtherDerived>                                   
      DipoleParticle& operator=(const Eigen::MatrixBase<OtherDerived> &other) {      
        Tvec::operator=(other);                                       
        return *this;                                                 
      }  
};

// Insert energy functions here!
// Major goal:
//   function: pairpot(particle &a, particle &b, point &vdist)

int main() {
  InputMap mcp("nemo.conf");
  Geometry::Sphere geo(mcp);
  Potential::Coulomb coulomb(mcp);

  cout << geo.info() << coulomb.info(25);

  DipoleParticle a, b;

  Point vdist = geo.vdist(a,b);
  double u = coulomb(a,b, geo.dist(a,b) );

  a.ranunit( slp_global ); // random unit vector
  double r = slp_global(); // random number [0:1[

  // Examples
  cout << "Squared distance = " << geo.sqdist(a,b) << "\n"; 
  cout << "Bjerrum          = " << coulomb.bjerrumLength() << "\n";
  cout << "Avogadro's number= " << pc::Nav << "\n";  // se "faunus/physconst.h"
}
