#ifndef FAU_POT_DEBYEHUCKEL_H
#define FAU_POT_DEBYEHUCKEL_H
#include "faunus/potentials/base.h"
namespace Faunus {
  /*! \brief Debye-Huckel potential
   *  \author Mikael Lund
   */
  class pot_debyehuckel : public pot_lj {
    private:
      double k;
    public:
      //! \param pot.lB Bjerrum length
      //! \param pot.eps L-J parameter
      //! \param pot.kappa Inverse Debye screening length
      pot_debyehuckel( pot_setup &pot ) : pot_lj(pot.eps/pot.lB) {
        f=pot.lB; 
        k=pot.kappa; 
        name+="/Debye-Huckel";
      };
      string info();
      //! \f$ \beta u/f = \frac{z_1z_2}{r}\exp(-\kappa r) + u_{lj}/f \f$
      //! \return Energy in kT/f (f=lB)
      inline double pairpot( const particle &p1, const particle &p2 ) {
        double r2=p1.sqdist(p2),
               r=sqrt(r2);
        return lj(p1,p2,r2) + p1.charge*p2.charge/r*exp(-k*r);
      }
  };
  string pot_debyehuckel::info() {
    std::ostringstream o;
    o << pot_lj::info()
      << "#   Bjerrum length    = " << f     << endl
      << "#   Debye length      = " << 1./k  << endl;
    return o.str();
  }
  typedef Faunus::pot_debyehuckel T_pairpot;
}
#endif
