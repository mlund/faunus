#ifndef FAU_POT_DEBYEHUCKELP3_H
#define FAU_POT_DEBYEHUCKELP3_H
#include "faunus/potentials/base.h"
namespace Faunus {
  /*! \brief Debye-Huckel potential for periodic boundry 
   *         conditions in 3D, it is extended to preform 
   *         under conditions of constant pressure.
   *         See class isobaric->markovmove.h
   *  \author Mikael Lund/Bjoern Persson
   *  \date Lund/Prag 2008
   *  \todo Remove pot_setop dependency
   */
  class pot_debyehuckelP3 : public pot_lj {
    private:
      double I,k;
    public:
      double box, halfbox;
      pot_debyehuckelP3( inputfile &in ) : pot_lj(in) {
        f=in.getflt("bjerrum",7.1);
        k=1/in.getflt("debyelen",1.1e6);
        if ( 1/k>=1e6) {
          I=in.getflt("ionicstr",0);
          k=sqrt( 4*std::acos(-1)*f*6.022e23/1e27*2*I );
        }
        box=in.getflt("boxlen");
        halfbox=box/2;
        eps=eps/f;
        name+="/Debye-Huckel w. minimum image";
      }

      //! \f$ \beta u/f = \frac{z_1z_2}{r}\exp(-\kappa r) + u_{lj}/f \f$
      //! \return Energy in units of kT/lB
      inline double pairpot( const particle &p1, const particle &p2 ) const {
        double r2=p1.sqdist(p2,box,halfbox), r=sqrt(r2);
        return lj(p1,p2,r2) + p1.charge*p2.charge/r*exp(-k*r);
      }

      inline double sqdist(const point &p1, const point &p2) {
        return p1.sqdist(p2,box,halfbox);
      }

      void setvolume(double vol) {
        box=pow(vol, 1./3);
        halfbox=box/2.;
      }

      string info() {
        std::ostringstream o;
        o << pot_lj::info()
          << "#   Bjerrum length    = " << f     << endl
          << "#   Debye length      = " << 1./k  << endl
          << "#   Ionic strength (M)= " << k*k*1e27/(8*std::acos(-1)*f*6.022e23) << endl;
        return o.str();
      }
  };
  class pot_debyehuckelP3trunk : public pot_lj_trunk {
    private:
      double I,k;
    public:
      double box, halfbox;
      pot_debyehuckelP3trunk( inputfile &in ) : pot_lj_trunk(in) {
        f=in.getflt("bjerrum",7.1);
        k=1/in.getflt("debyelen",1.1e6);
        if ( 1/k>=1e6) {
          I=in.getflt("ionicstr",0);
          k=sqrt( 4*std::acos(-1)*f*6.022e23/1e27*2*I );
        }
        box=in.getflt("boxlen");
        halfbox=box/2;
        eps=eps/f;
        name+="/Debye-Huckel w. minimum image, augmented with truncated LJ";
      }

      //! \f$ \beta u/f = \frac{z_1z_2}{r}\exp(-\kappa r) + u_{lj}/f \f$
      //! \return Energy in units of kT/lB
      inline double pairpot( const particle &p1, const particle &p2 ) const {
        double r2=p1.sqdist(p2,box,halfbox), r=sqrt(r2);
        return lj(p1,p2,r2) + p1.charge*p2.charge/r*exp(-k*r);
      }

      inline double sqdist(const point &p1, const point &p2) {
        return p1.sqdist(p2,box,halfbox);
      }

      void setvolume(double vol) {
        box=pow(vol, 1./3);
        halfbox=box/2.;
      }

      string info() {
        std::ostringstream o;
        o << pot_lj_trunk::info()
          << "#   Bjerrum length    = " << f     << endl
          << "#   Debye length      = " << 1./k  << endl
          << "#   Ionic strength (M)= " << k*k*1e27/(8*std::acos(-1)*f*6.022e23) << endl;
        return o.str();
      }
  };
}
#endif
