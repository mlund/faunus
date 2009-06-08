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
    public:
      string name;
      double f;
      pot_hsminimage( inputfile &in ) {
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
        double r=p1.dist(p2,box,halfbox), u=p1.charge*p2.charge/r;
        return (r<p1.radius+p2.radius) ? u+200. : u;
      }
      inline double sqdist(const point &p1, const point &p2) {
        return p1.sqdist(p2,box,halfbox);
      }
      string info() {
        std::ostringstream o;
        o << "#   Type              = " << name << std::endl
          << "#   Bjerrum length    = " << f << std::endl
          << "#   Image length      = " << box << std::endl;
        return o.str();
      }
  };

  /*!
   * \brief Coulomb/Hardsphere potential with minimum image in XY, hard in Z.
   * \author Mikael Lund
   * \date 2009
   */
  class pot_hsminimageXY {
    private:
      double halfbox,box;
    public:
      string name;
      double f;
      pot_hsminimageXY( inputfile &in ) {
        name+="Coulomb/Hardsphere w. minimum image (XY, only)";
        f=in.getflt("bjerrum",7.1);
        box=in.getflt("boxlen");
        halfbox=box/2.;
      }
      void setvolume(double vol) {
        box=pow(vol, 1./3);
        halfbox=box/2.;
      }
      inline double pairpot(const particle &p1, const particle &p2) {
        double r=sqrt(sqdist(p1,p2)), u=p1.charge*p2.charge/r;
        return (r<p1.radius+p2.radius) ? u+200. : u;
      }
      inline double sqdist(const particle &p1, const particle &p2) {
        double dz=p1.z-p2.z;
        double dx=std::abs(p1.x-p2.x);
        double dy=std::abs(p1.y-p2.y);
        if (dx>halfbox) dx-=box;
        if (dy>halfbox) dy-=box;
        return dx*dx + dy*dy + dz*dz;
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
