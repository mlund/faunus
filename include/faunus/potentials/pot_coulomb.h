#ifndef FAU_POT_COULOMB_H
#define FAU_POT_COULOMB_H
#include "faunus/potentials/base.h"
namespace Faunus {
  /*! \brief Coulomb potential
   *  \author Mikael Lund
   *  \date Prague, 2007
   *
   *  This potential class calculates the solvent screened pair
   *  potential between two charges. The Faunus::inputfile
   *  constructor argument looks for the Bjerrum length, "bjerrum".
   *  If not found that of water at 298K is used (7.1 Angstrom).
   *
   *  Faunus::inputfile keywords:
   *  \li bjerrum - The Bjerrum length, \f$l_B\f$ [AAngstrom]
   *  \li LJeps - Lennard-Jones interaction parameter [kT] (see Faunus::pot_lj::eps)
   */
  class pot_coulomb : public pot_lj {
    public:
      pot_coulomb(const inputfile &in) : pot_lj(in) {
        f=in.getflt("bjerrum",7.1);
        eps=eps/f;
        name+="/Coulomb";
      }
      /*! \brief Return Coulomb energy between a pair of particles
       *  \return Energy in units of \f$kT/l_B\f$ (lB=f). \f[ \beta u/l_B = \frac{z_1 z_2}{r} + \frac{u_{lj}}{l_B} \f]
       */
      inline double pairpot(const particle &p1, const particle &p2) {
        register double r2=p1.sqdist(p2);
        return lj(p1,p2,r2) + p1.charge*p2.charge/sqrt(r2);
      }

      string info() {
        std::ostringstream o;
        o << pot_lj::info()
          << "#   Bjerrum length    = " << f << endl;
        return o.str();
      }
  };
}
#endif
