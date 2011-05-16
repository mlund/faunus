#ifndef FAU_POT_SILKSPIDER_H
#define FAU_POT_SILKSPIDER_H
#include "faunus/potentials/base.h"
namespace Faunus {
  class pot_silkspider : public pot_harmonic {
    private:
      double halfbox, box, I, sticky_r, sticky_u;
    public:
      double f, kappa;
      string name;

      pot_silkspider( inputfile &in ) : pot_harmonic(in) {
        name="r12 + Debye-Huckel w. minimum image and hydrophobic stickiness";
        f=in.getflt("bjerrum",7.1);
        box=in.getflt("boxlen");
        halfbox=box/2;
        I=in.getflt("ionicstr",0);
        sticky_r=in.getflt("sticky_r", 0);   // sticky threshold between *surfaces*
        sticky_u=in.getflt("sticky_u", 0)/f; // sticky energy in units of kT
        k=sqrt( 8*pyc.pi*f*pyc.Nav/1e27*I );
      }

      void setvolume(double vol) {
        box=pow(vol, 1./3);;
        halfbox=box/2;
      }

      inline double pairpot(const particle &p1, const particle &p2) {
        double r2=p1.sqdist_mi_xyz(p2,box,halfbox), s=p1.radius+p2.radius, a=s*s/r2;
        s=a*a*a;
        double u=s*s/f;
        r2=sqrt(r2);
        if (p1.hydrophobic==true)
          if (p2.hydrophobic==true)
            if (r2-p1.radius-p2.radius < sticky_r)
              u+=sticky_u;
        return u + p1.charge*p2.charge*exp(-k*r2)/r2;
      }

      inline double sqdist(const point &p1, const point &p2) {
        return p1.sqdist_mi_xyz(p2,box,halfbox);
      }

      inline double bond(particle &p1, particle &p2) {
        req=p1.radius+p2.radius+2.;
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
          << "#   Bjerrum length    = " << f << " AA\n"
          << "#   Debye length      = " << 1/k << " AA\n"
          << "#   Ionic strength (M)= " << I << " mol/l\n"
          << "#   Sticky threshold  = " << sticky_r << " AA\n"
          << "#   Sticky energy     = " << sticky_u*f << " kT\n"
          << "#   Image length      = " << box << " AA" << endl;
        return o.str();
      }
  };
} //namespace
#endif
