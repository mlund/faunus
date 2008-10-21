#ifndef FAU_POT_HSMINIMAGE_H
#define FAU_POT_HSMINIMAGE_H
#include "faunus/potentials/base.h"
namespace Faunus {
  /*!
   * \brief Coulomb/Hardsphere potential with minimum image.
   * \author Mikael Lund
   * \date 2008
   */
  class pot_hsminimage {
    private:
      double halfbox,box;
      string name;
    public:
      double f;
      pot_hsminimage( const inputfile &in ) {
        name+="Coulomb/Hardsphere w. minimum image";
        f=in.getflt("bjerrum",7.1);
        box=in.getflt("boxlen");
        halfbox=box/2.;
      }
      void setvolume(double vol) {
        box=pow(vol, 1./3);
        halfbox=box/2.;
      }
      inline double pairpot(const particle &p1, const particle &p2) {
        double r2=p1.sqdist(p2,box,halfbox), s=p1.radius+p2.radius;
        return (r2<s*s) ? 20000. : p1.charge*p2.charge/sqrt(r2);
      }
      string info() {
        std::ostringstream o;
        o << "#   Type              = " << name << std::endl
          << "#   Bjerrum length    = " << f << std::endl
          << "#   Image length      = " << box << std::endl;
        return o.str();
      }
  };
}
#endif
