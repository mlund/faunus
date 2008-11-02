#ifndef FAU_POT_BARECOULOMB_H
#define FAU_POT_BARECOULOMB_H
#include "faunus/potentials/base.h"
namespace Faunus {
  /*! \brief Naked coulomb potential
   *  \author Mikael Lund
   *  \date Prague, 2008
   *  \note No vdW or hardsphere added.
   */
  class pot_barecoulomb {
    public:
      double f;
      string name;
      pot_barecoulomb (inputfile &in) {
        f=in.getflt("bjerrum",7.1);
        name="Bare coulomb potential";
      }
      /*! Return Coulomb energy between a pair of particles
       *
       *  \return Energy in units of kT/f (f=lB).
       *  \f[ \beta u/l_B = \frac{z_1 z_2}{r} \f]
       */
      inline double pairpot(const particle &p1, const particle &p2) {
        return p1.charge*p2.charge/p1.dist(p2);
      }
      virtual string info() {
        std::ostringstream o;
        o << "#   Type              = " << name << std::endl
          << "#   Bjerrum length    = " << f << std::endl;
        return o.str();
      }
  };
}
#endif
