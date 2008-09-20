#ifndef FAU_POT_HSCOULOMB_H
#define FAU_POT_HSCOULOMB_H
#include "faunus/potentials/base.h"
namespace Faunus {
  /*! \brief Coulomb potential
   *  \author Mikael Lund
   *  \date Prague, 2007
   */
  class pot_hscoulomb : public pot_hs {
    public:
      pot_hscoulomb (const inputfile &in) : pot_hs() {
        f=in.getflt("bjerrum",7.1);
        name+="/Coulomb";
      }
      /*! \brief Return Coulomb energy between a pair of particles
       *  \return Energy in units of kT/f (f=lB).
       *  \f$ \beta u/f = \frac{z_1 z_2}{r} + u_{HS}/f \f$
       */
      inline double pairpot(const particle &p1, const particle &p2) {
        register double r2=p1.sqdist(p2);
        return hs(p1,p2,r2) + p1.charge*p2.charge/sqrt(r2);
      }
      string info() {
        std::ostringstream o;
        o << pot_hs::info()
          << "#   Bjerrum length    = " << f << endl;
        return o.str();
      }
  };
}
#endif
