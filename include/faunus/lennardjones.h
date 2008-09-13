#ifndef FAU_LENNARDJONES_H
#define FAU_LENNARDJONES_H

#include <sstream>
#include "faunus/point.h"

namespace Faunus {
  /*!
   *  \brief Lennard-Jones potential
   *  \author Mikael Lund
   *  \date Prague, 2007
   */
  class pot_lj {
    public:
      string name; //!< Arbitrary name
      string cite; //!< Litterature reference
      double eps;  //!< 4*Lennard-Jones interaction parameter (kT)
      double f;    //!< Factor to convert to kT (used after energy summations)
      pot_lj(double epsilon) {
        eps=epsilon;
        name="LJ12-6";
        f=1;
      }
      /*!
       *  L-J pair energy.
       *  \f$ u_{lj} = \epsilon \left ( \frac{\sigma}{r^{12}} - \frac{\sigma}{r^6} \right ) \f$
       *  \param r2 Squared distance between particle 1 and 2.
       */
      inline double lj(const particle &p1, const particle &p2, double &r2) {
        register double x=p1.radius+p2.radius,u=x*x/r2;
        x=u*u*u;
        return (x*x-x)*eps;
      }
      inline void lj(const particle &p1, const particle &p2, const double &r2, double &u) {
        register double s=p1.radius+p2.radius, a=s*s/r2;
        s=a*a*a;
        u+=(s*s-s)*eps;
      }
      virtual void setvolume(double) {}; //!< Function to specify volume for fluctuating periodic boundaries
      virtual string info() {
        std::ostringstream o;
        o << "#   Type              = " << name << std::endl
          << "#   4*LJ epsilon (kT) = " << eps*f << std::endl;
        if (cite.length()!=0)
          o << "#   Reference         = " << cite << std::endl;
        return o.str();
      }
  };

  /*!
   *  \brief Hardsphere potential
   *  \author Mikael Lund
   *  \date Prague, 2008
   */
  class pot_hs {
    public:
      string name; //!< Arbitrary name
      string cite; //!< Litterature reference
      double f;    //!< Factor to convert to kT (used after energy summations)
      pot_hs(double dummy=0) {
        name="Hardsphere";
        f=1;
      }
      inline double hs(const particle &p1, const particle &p2, const double &r2) {
        register double s=p1.radius+p2.radius;
        return (r2>s*s) ? 0 : 99999.;
      }
      inline void hs(const particle &p1, const particle &p2, const double &r2, double &u) { u+=hs(p1,p2,r2); }
      virtual void setvolume(double) {}; //!< Function to specify volume for fluctuating periodic boundaries
      virtual string info() {
        std::ostringstream o;
        o << "#   Type              = " << name << std::endl;
        if (cite.length()!=0)
          o << "#   Reference         = " << cite << std::endl;
        return o.str();
      }
  };
}//namespace
#endif
