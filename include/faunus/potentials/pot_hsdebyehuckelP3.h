#ifndef FAU_POT_HSDEBYEHUCKELP3_H
#define FAU_POT_HSDEBYEHUCKELP3_H
#include "faunus/potentials/base.h"
#include "faunus/fortran.h"

namespace Faunus {
  /*! \brief Debye-Huckel potential for periodic boundry 
   *         conditions in 3D with hard sphere repulsion
   *  \author Anil Kurut
   *  \date Lund 2010
   */
  class pot_hsdebyehuckelP3 {
    protected:
      double I,k;
    public:
      double box, halfbox, f;
      string name;
      pot_hsdebyehuckelP3( inputfile &in ) {
        f=in.getflt("bjerrum",7.1);
        k=1/in.getflt("debyelen",1.1e6);
        if ( 1/k>=1e6) {
          I=in.getflt("ionicstr",0);
          k=sqrt( 8*pyc.pi*f*pyc.Nav/1e27*I );
        }
        box=in.getflt("boxlen");
        halfbox=box/2;
        name+="/Debye-Huckel w. hard sphere and minimum image";
      }

      double inline f_exp(double x) {
        static double e=2.718281828;
        int tmp = (*(1 + (int *) &e));
        int tmp2 = (int)(x * (tmp - 1072632447) + 1072632447);
        double p = 0.0;
        *(1 + (int * ) &p) = tmp2;
        return p;
      }

      //! \f$ \beta u/f = \frac{z_1z_2}{r(1+\kappa a)}\exp(-\kappa (r-a) ) + u_{hs} \f$
      //! \return Energy in units of kT/lB
      inline double pairpot( const particle &p1, const particle &p2 ) {
        double r=p1.sqdist(p2,box,halfbox),
               a=p1.radius+p2.radius;
        if (r<a*a)
          return 5000.;
        r=sqrt(r);
        return p1.charge*p2.charge / (r*(1.+k*a)) * f_exp( k*(a-r) );
      }

      inline double sqdist(const point &p1, const point &p2) {
        return p1.sqdist(p2,box,halfbox);
      }

      void setvolume(double vol) {
        box=pow(vol, 1./3.);
        halfbox=0.5*box;
      }

      int anint(double x) const {
        return int(x>0. ? x+.5 : x-.5);
      }

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
          << "#   Current boxlength = " << box << endl
          << "#   Ionic strength (M)= " << k*k*1e27/(8*pyc.pi*f*pyc.Nav) << endl;
        return o.str();
      }
  };
} // namespace
#endif
