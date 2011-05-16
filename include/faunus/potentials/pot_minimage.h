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
        double r2=p1.sqdist_mi_xyz(p2,box,halfbox);
        return lj(p1,p2,r2) + p1.charge*p2.charge/sqrt(r2);
      }
      inline double sqdist(const point &p1, const point &p2) {
        return p1.sqdist_mi_xyz(p2,box,halfbox);
      }
      int anint(double x) const { return int(x>0. ? x+.5 : x-.5); }
      void boundary(point &a) const {
        if (std::abs(a.x)>halfbox) a.x-=box*anint(a.x/box);
        if (std::abs(a.y)>halfbox) a.y-=box*anint(a.y/box);
        if (std::abs(a.z)>halfbox) a.z-=box*anint(a.z/box);
      }

      string info() {
        std::ostringstream o;
        o << pot_lj::info()
          << "#   Bjerrum length    = " << f << endl
          << "#   Image length      = " << box << endl;
        return o.str();
      }
  };

  class pot_r12minimageXY : public pot_harmonic {
    private:
      double halfbox,box;
    public:
      double f;       //!< Factor to convert to kT
      string name;

      pot_r12minimageXY( inputfile &in ) : pot_harmonic(in) {
        name="r12 + Coulomb w. minimum image in XY-directions";
        f=in.getflt("bjerrum",7.1);
        box=in.getflt("boxlen");
        halfbox=box/2;
      }

      void setvolume(double vol) {
        box=pow(vol, 1./3);
        halfbox=box/2;
      }

      inline double pairpot(const particle &p1, const particle &p2) {
        double r2=sqdist(p1,p2), s=p1.radius+p2.radius, a=s*s/r2;
        s=a*a*a;
        return s*s/f + p1.charge*p2.charge/sqrt(r2);
      }

      inline double sqdist(const point &p1, const point &p2) {
        double dz=p1.z-p2.z;
        double dx=std::abs(p1.x-p2.x);
        double dy=std::abs(p1.y-p2.y);
        if (dx>halfbox) dx-=box;
        if (dy>halfbox) dy-=box;
        return dx*dx + dy*dy + dz*dz;
      }

      /*!
       * Calculated the bonded + non-bonded energy between two particles
       * connected with a harmonic bond. To avoid numerical trouble with
       * the non-electrostatic, non-bonded interactions, the radii of
       * the two particles are temporarily reduced.
       */
      inline double bond(particle &p1, particle &p2) {
        double r=sqrt(sqdist(p1,p2));
        double u=harmonicbond(p1,p2,r)/f;
        p1.radius*=0.5;
        p2.radius*=0.5;
        u+=pairpot(p1,p2);
        p1.radius*=2;
        p2.radius*=2;
        return u;
      }

      string info() {
        std::ostringstream o;
        o << "#     Name              = " << name << endl
          << "#     Bjerrum length    = " << f << endl
          << "#     Image length (XY) = " << box << endl
          << pot_harmonic::info();
        return o.str();
      }
  };

  class pot_r12minimage : public pot_harmonic {
    private:
      double halfbox,box;
    public:
      double f;
      string name;
      pot_r12minimage( inputfile &in ) : pot_harmonic(in) {
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
        double r2=p1.sqdist_mi_xyz(p2,box,halfbox), s=p1.radius+p2.radius, a=s*s/r2;
        s=a*a*a;
        return s*s/f + p1.charge*p2.charge/sqrt(r2);
      }
      inline double sqdist(const point &p1, const point &p2) {
        return p1.sqdist_mi_xyz(p2,box,halfbox);
      }
      inline double bond(particle &p1, particle &p2) {
        double r=sqrt(sqdist(p1,p2));
        double u=harmonicbond(p1,p2,r)/f;
        p1.radius*=0.5;
        p2.radius*=0.5;
        u+=pairpot(p1,p2);
        p1.radius*=2;
        p2.radius*=2;
        return u;
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
