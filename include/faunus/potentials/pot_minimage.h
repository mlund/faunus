#ifndef FAU_POT_MINIMAGE_H
#define FAU_POT_MINIMAGE_H
#include "faunus/potentials/base.h"
namespace Faunus {
  /*!
   * \brief Coulomb pot. with minimum image.
   * \author Mikael Lund
   * \date 2007
   */
  class pot_minimage : public pot_lj {
    private:
      double invbox,box;
    public:
      pot_minimage(pot_setup &pot) : pot_lj(pot.eps/pot.lB) {
        f=pot.lB;
        box=pot.box;
        invbox=1./box;
        name+="/Coulomb w. minimum image";
      }
      void setvolume(double vol) {
        box=pow(vol, 1./3);;
        invbox=1./box;
      }
      inline double pairpot(const particle &p1, const particle &p2) {
        register double r2=p1.sqdist(p2,box,invbox);
        return lj(p1,p2,r2) + p1.charge*p2.charge/sqrt(r2);
      }
      string info() {
        std::ostringstream o;
        o << pot_lj::info()
          << "#   Bjerrum length    = " << f << endl
          << "#   Image length      = " << box << endl;
        return o.str();
      }
  };
}
#endif
