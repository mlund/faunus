#ifndef FAU_POT_BASE_H
#define FAU_POT_BASE_H
#include "faunus/common.h"
#include "faunus/point.h"
#include "faunus/inputfile.h"

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
      double eps;  //!< Lennard-Jones interaction parameter (kT)
      double f;    //!< Factor to convert to kT (used after energy summations)
      pot_lj(const inputfile &in) {
        eps=in.getflt("LJeps", 2);
        name="LJ12-6";
        f=1;
      }
      pot_lj(double epsilon) {
        eps=epsilon;
        name="LJ12-6";
        f=1;
      }
      /*!
       *  \param r2 Squared distance between particle p1 and p2.
       *  \returns Interaction energy in units of kT,
       *           \f[ \beta u_{lj} = \epsilon_{lj} \left ( \frac{\sigma}{r^{12}} - \frac{\sigma}{r^6} \right ) \f]
       *           \f[ \sigma = \frac{\sigma_1+\sigma_2}{2}\f]
       */
      inline double lj(const particle &p1, const particle &p2, double &r2) const {
        register double x=p1.radius+p2.radius,u=x*x/r2;
        x=u*u*u;
        return (x*x-x)*eps;
      }
      inline void lj(const particle &p1, const particle &p2, const double &r2, double &u) const {
        register double s=p1.radius+p2.radius, a=s*s/r2;
        s=a*a*a;
        u+=(s*s-s)*eps;
      }
      virtual void setvolume(double) {}; //!< specify volume for fluctuating periodic boundaries
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
      inline double hs(const particle &p1, const particle &p2, const double &r2) const {
        register double s=p1.radius+p2.radius;
        return (r2>s*s) ? 0 : 99999.;
      }
      inline void hs(const particle &p1, const particle &p2, const double &r2, double &u) const {
        u+=hs(p1,p2,r2);
      }
      virtual void setvolume(double) {}; //!< specify volume for fluctuating periodic boundaries
      virtual string info() {
        std::ostringstream o;
        o << "#   Type              = " << name << std::endl;
        if (cite.length()!=0)
          o << "#   Reference         = " << cite << std::endl;
        return o.str();
      }
  };
}
#endif
