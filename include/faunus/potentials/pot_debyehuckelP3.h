#ifndef FAU_POT_DEBYEHUCKELP3_H
#define FAU_POT_DEBYEHUCKELP3_H
#include "faunus/potentials/base.h"
namespace Faunus {
  /*! \brief Debye-Huckel potential for periodic boundry 
   *         conditions in 3D, it is extended to preform 
   *         under conditions of constant pressure.
   See class isobaric->markovmove.h
   *  \author Mikael Lund/Bjoern Persson
   *  \date Lund/Prag 2008
   */
  class pot_debyehuckelP3 : public pot_lj {
    private:
      double k;
    public:
      double box, invbox;
      //! \param pot.lB Bjerrum length
      //! \param pot.eps L-J parameter
      //! \param pot.kappa Inverse Debye screening length
      pot_debyehuckelP3( pot_setup &pot ) : pot_lj(pot.eps/pot.lB) {
        f=pot.lB; 
        k=pot.kappa;
        box=pot.box;
        invbox=1./box;
        name+="/Debye-Huckel w. minimum image";
      };
      string info();
      //! \f$ \beta u/f = \frac{z_1z_2}{r}\exp(-\kappa r) + u_{lj}/f \f$
      //! \return Energy in kT/f (f=lB)
      inline double pairpot( const particle &p1, const particle &p2 ) {
        register double r2=p1.sqdist(p2,box,invbox), r=sqrt(r2);
        return lj(p1,p2,r2) + p1.charge*p2.charge/r*exp(-k*r);
      }
      void setvolume(double vol) {
        box=pow(vol, 1./3);;
        invbox=1./box;
      }
  }; 
  string pot_debyehuckelP3::info() {
    std::ostringstream o;
    o << pot_lj::info()
      << "#   Bjerrum length    = " << f     << endl
      << "#   Debye length      = " << 1./k  << endl;
    return o.str();
  }
  typedef Faunus::pot_debyehuckelP3 T_pairpot;
}
#endif
