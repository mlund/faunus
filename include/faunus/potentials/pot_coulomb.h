#ifndef FAU_POT_COULOMB_H
#define FAU_POT_COULOMB_H
#include "faunus/potentials/base.h"
namespace Faunus {
  /*! \brief Coulomb potential
   *  \author Mikael Lund
   *  \date Prague, 2007
   *  \todo Remove the hackish typedef to allow pot_netz to enherit
   */
  class pot_coulomb : public pot_lj {
    public:
      /*! \param pot.lB Bjerrum length
       *  \param pot.eps L-J epsilon parameter (in kT) */
      pot_coulomb ( pot_setup &pot) : pot_lj(pot.eps/pot.lB) {
        f=pot.lB;
        name+="/Coulomb";
      }
      /*! \brief Return Coulomb energy between a pair of particles
       *  \return Energy in units of kT/f (f=lB).
       *  \f$ \beta u/f = \frac{z_1 z_2}{r} + u_{LJ}/f \f$
       */
      inline double pairpot(const particle &p1, const particle &p2) {
        register double r2=p1.sqdist(p2);
        return lj(p1,p2,r2) + p1.charge*p2.charge/sqrt(r2);
      }
      string info();
  };
  string pot_coulomb::info() {
    std::ostringstream o;
    o << pot_lj::info()
      << "#   Bjerrum length    = " << f << endl;
    return o.str();
  }
#ifndef _POT_NETZ_H
  typedef Faunus::pot_coulomb T_pairpot;
#endif
}
#endif
