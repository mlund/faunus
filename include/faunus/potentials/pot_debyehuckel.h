#ifndef FAU_POT_DEBYEHUCKEL_H
#define FAU_POT_DEBYEHUCKEL_H
#include "faunus/potentials/base.h"
namespace Faunus {
  /*! \brief Debye-Huckel potential
   *  \author Mikael Lund
   *  
   *  Faunus::inputfile is scanned for "bjerrum", "debyelen", "LJeps".
   */
  class pot_debyehuckel : public pot_lj {
    public:
      double I,k;
      pot_debyehuckel( inputfile &in ) : pot_lj(in) {
        name+="/Debye-Huckel";
        f=in.getflt("bjerrum",7.1);
        k=1./in.getflt("debyelen",1.1e4);
        if ( 1/k>=1e4) {
          I=in.getflt("ionicstr",0);
          k=sqrt( 4*std::acos(-1)*f*6.022e23/1e27*2*I );
        }
        eps=eps/f;
      }
      //! \f$ \beta u/f = \frac{z_1z_2}{r}\exp(-\kappa r) + u_{lj}/f \f$
      //! \return Energy in kT/f (f=lB)
      inline double pairpot( const particle &p1, const particle &p2 ) {
        double r2=p1.sqdist(p2),
               r=sqrt(r2);
        return lj(p1,p2,r2) + p1.charge*p2.charge/r*exp(-k*r);
      }
      inline double sqdist(const point &p1, const point &p2) {
        return p1.sqdist(p2);
      }
      string info() {
        std::ostringstream o;
        o << pot_lj::info()
          << "#   Bjerrum length    = " << f << endl
          << "#   Kappa             = " << k << endl
          << "#   Debye length      = " << 1./k << endl
          << "#   Ionic strength (M)= " << k*k*1e27/(8*std::acos(-1)*f*6.022e23) << endl;
        return o.str();
      }
  };
}
#endif
