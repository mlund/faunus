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
      double halfbox,box;
    public:
      pot_minimage( inputfile &in ) : pot_lj(in) {
        name+="/Coulomb w. minimum image";
        f=in.getflt("bjerrum",7.1);
        box=in.getflt("boxlen");
        halfbox=box/2;
        eps=eps/f;
      }
      void setvolume(double vol) {
        box=pow(vol, 1./3);;
        halfbox=box/2;
      }
      inline double pairpot(const particle &p1, const particle &p2) {
        double r2=p1.sqdist(p2,box,halfbox);
        return lj(p1,p2,r2) + p1.charge*p2.charge/sqrt(r2);
      }
      inline double sqdist(const point &p1, const point &p2) {
        return p1.sqdist(p2,box,halfbox);
      }
      string info() {
        std::ostringstream o;
        o << pot_lj::info()
          << "#   Bjerrum length    = " << f << endl
          << "#   Image length      = " << box << endl;
        return o.str();
      }
  };

  class pot_r12minimage {
    private:
      double halfbox,box;
    public:
      double f;
      string name;
      pot_r12minimage( inputfile &in ) {
        name="r12 + Coulomb w. minimum image";
        f=in.getflt("bjerrum",7.1);
        box=in.getflt("boxlen");
        halfbox=box/2;
      }
      void setvolume(double vol) {
        box=pow(vol, 1./3);;
        halfbox=box/2;
      }
      inline double pairpot(const particle &p1, const particle &p2) {
        double r2=p1.sqdist(p2,box,halfbox), s=p1.radius+p2.radius, a=s*s/r2;
        s=a*a*a;
        return s*s/f + p1.charge*p2.charge/sqrt(r2);
      }
      inline double sqdist(const point &p1, const point &p2) {
        return p1.sqdist(p2,box,halfbox);
      }
      string info() {
        std::ostringstream o;
        o << "#   Name              = " << name << endl
          << "#   Bjerrum length    = " << f << endl
          << "#   Image length      = " << box << endl;
        return o.str();
      }
  };
}
#endif
