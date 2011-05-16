#ifndef FAU_POT_HSHAMAKERDH_H
#define FAU_POT_HSHAMAKERDH_H
#include "faunus/potentials/base.h"
#include "faunus/fortran.h"

namespace Faunus {
  /*! \brief Debye-Huckel potential with hard sphere repulsion and hamaker attraction between spheres
   *         for periodic boundry conditions in 3D, 
   *         it is extended to preform 
   *         under conditions of constant pressure.
   *         
   *  \author Mikael Lund/Bjoern Persson
   *  \date Lund/Prag 2008
   *  \todo Remove pot_setop dependency
   */
  class pot_hsHamakerDH {
    protected:
      double I,k;
    public:
      string name;
      double box, halfbox, f, A;
      pot_hsHamakerDH( inputfile &in ) {
        f=in.getflt("bjerrum",7.1);
        k=1/in.getflt("debyelen",1.1e6);
        A=in.getflt("hamaker",9)/f;      // in units of kT/lB
        if ( 1/k>=1e6) {
          I=in.getflt("ionicstr",0);
          k=sqrt( 4*std::acos(-1.)*f*6.022e23/1e27*2*I );
        }
        box=in.getflt("boxlen");
        halfbox=box/2;
        name+="/Debye-Huckel w. hard spehere, hamaker, and minimum image";
      }

      double inline f_exp(double x) {
        static double e=2.718281828;
        int tmp = (*(1 + (int *) &e));
        int tmp2 = (int)(x * (tmp - 1072632447) + 1072632447);
        double p = 0.0;
        *(1 + (int * ) &p) = tmp2;
        return p;
      }

      //! \f$ \beta u/f = \frac{z_1z_2}{r}\exp(-\kappa r) \f$
      //! \return Energy in units of kT/lB
      inline double pairpot( const particle &p1, const particle &p2 ) {
        double r2=p1.sqdist_mi_xyz(p2,box,halfbox),
               r=sqrt(r2), 
               x=p1.radius+p2.radius,
               //x2=x*x,
               y=p1.radius*p2.radius,
               //y2=y*y,
               //z=p1.radius-p2.radius,
               //z2=z*z,
               u=p1.charge * p2.charge / r * f_exp(-k*r) - A*y/(6*(r-x+2)*x);
        //u=p1.charge*p2.charge/r*f_exp(-k*r)-A/6*((2*y/(r2-x2+2))+(2*y/(r2-z2))+std::log((r2-x2+2)/(r2-z2)));

        //if (r>x)
          //std::cout << "u = " << u << std::endl;
    
        return (r<x) ? u+2000 : u;
      }

      inline double sqdist(const point &p1, const point &p2) {
        return p1.sqdist_mi_xyz(p2,box,halfbox);
      }

      void setvolume(double vol) {
        box=pow(vol, 1./3.);
        halfbox=box/2.;
      }

      int anint(double x) const { return int(x>0. ? x+.5 : x-.5); }
      void boundary(point &a) const {
        if (std::abs(a.x)>halfbox) a.x-=box*anint(a.x/box);
        if (std::abs(a.y)>halfbox) a.y-=box*anint(a.y/box);
        if (std::abs(a.z)>halfbox) a.z-=box*anint(a.z/box);
      }

      string info() {
        std::ostringstream o;
        o << "#   Potential type    = " << name << endl
          << "#   Bjerrum length    = " << f     << endl
          << "#   Debye length      = " << 1./k  << endl
          << "#   Ionic strength (M)= " << k*k*1e27/(8*std::acos(-1.)*f*6.022e23) << endl
          << "#   Hamaker constant  = " << A*f << endl;
        return o.str();
      }
  };
}
#endif
