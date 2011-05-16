#ifndef FAU_POT_BASE_H
#define FAU_POT_BASE_H
#include "faunus/common.h"
#include "faunus/point.h"
#include "faunus/inputfile.h"
#include "faunus/physconst.h"

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
      pot_lj(inputfile &in) {
        eps=in.getflt("LJeps", -1 );          //!< DEPRECATED!! Use "lj_epsilon" instead
        if ( eps == -1 ) {
        	/*!< "lj_epsilon" is the new option for "eps" definition.
        	 *   The value is multiplied by 4 to follow the standard notation
        	 *   of the LJ potential
        	 */
        	eps=in.getflt("lj_epsilon", 0.2 );
        	eps*=4;
        }
        name="LJ12-6";
        f=1;
      }
      pot_lj(double epsilon) {
        eps=epsilon;
        name="LJ12-6";
        f=1;
      }
      /*!
       *  \param p1 First particle
       *  \param p2 Second particle
       *  \param r2 Squared distance between particle p1 and p2.
       *  \returns Interaction energy in units of kT,
       *           \f[ \beta u_{lj} = \epsilon_{lj} \left ( \frac{\sigma}{r^{12}} - \frac{\sigma}{r^6} \right ) \f]
       *           \f[ \sigma = \frac{\sigma_1+\sigma_2}{2}\f]
       */
      inline double lj(const particle &p1, const particle &p2, const double &r2) const {
        double x=p1.radius+p2.radius,u=x*x/r2;
        x=u*u*u;
        return (x*x-x)*eps;
      }
      inline void lj(const particle &p1, const particle &p2, const double &r2, double &u) const {
        double s=p1.radius+p2.radius, a=s*s/r2;
        s=a*a*a;
        u+=(s*s-s)*eps;
      }
      virtual void setvolume(double) {}  //!< specify volume for fluctuating periodic boundaries
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
        double s=p1.radius+p2.radius;
        return (r2>s*s) ? 0 : 9999.;
      }
      /*
      inline void hs(const particle &p1, const particle &p2, const double &r2, double &u) const {
        u+=hs(p1,p2,r2);
      }*/
      virtual void setvolume(double) {}  //!< specify volume for fluctuating periodic boundaries
      virtual string info() {
        std::ostringstream o;
        o << "#   Type              = " << name << std::endl;
        if (cite.length()!=0)
          o << "#   Reference         = " << cite << std::endl;
        return o.str();
      }
  };
  class pot_lj_trunk {
    public:
      string name; //!< Arbitrary name
      string cite; //!< Litterature reference
      double eps;  //!< Lennard-Jones interaction parameter (kT)
      double f;    //!< Factor to convert to kT (used after energy summations)
      double third2;
      pot_lj_trunk(inputfile &in) {
        eps=in.getflt("LJeps", 2);
        name="Truncated and shifted LJ12-6";
        f=1;
        third2=pow(2.,1./3.);
      }
      pot_lj_trunk(double epsilon) {
        eps=epsilon;
        name="Truncated and shifted LJ12-6";
        f=1;
        third2=pow(2.,1./3.);
      }
      /*!
       * \param p1 First particle
       * \param p2 Second particle
       * \param r2 Squared distance between particle p1 and p2.
       * \returns Interaction energy in units of kT,
       *           \f[ \beta u_{lj} = \epsilon_{lj} \left ( \frac{\sigma}{r^{12}} - \frac{\sigma}{r^6} \right ) \f]
       *           \f[ \sigma = \frac{\sigma_1+\sigma_2}{2}\f]
       */
      inline double lj(const particle &p1, const particle &p2, const double &r2) const {

        double x=p1.radius+p2.radius;
        if (r2>x*x*third2)
          return 0;
        else {
          double u=x*x/r2;
          x=u*u*u;
          return (x*x-x+0.25)*eps;
        }
      }
      inline void lj(const particle &p1, const particle &p2, const double &r2, double &u) const {
        double s=p1.radius+p2.radius, a=s*s/r2;
        s=a*a*a;
        u+=(s*s-s)*eps;
      }
      virtual void setvolume(double) {}  //!< specify volume for fluctuating periodic boundaries
      virtual string info() {
        std::ostringstream o;
        o << "#   Type              = " << name << std::endl
          << "#   4*LJ epsilon (kT) = " << eps*f << std::endl;
        if (cite.length()!=0)
          o << "#   Reference         = " << cite << std::endl;
        return o.str();
      }
  };

  class pot_harmonic {
    public:
      string name;    //!< Arbitrary name
      double k;       //!< Force constant for bonds (if any)
      double req;     //!< Equilibrium distance (in AA)
      pot_harmonic(inputfile &in) {
        k=in.getflt("harmonic_k", 0.3);
        req=in.getflt("harmonic_req", 0);
        name="Harmonic spring interaction";
      }

      double harmonicbond(double r) {
        double dr=r-req;
        return k*dr*dr;
      }

      double harmonicbond(particle &p1, particle &p2, double r) {
        double dr=r-req;
        return k*dr*dr;
      }

      string info() {
        std::ostringstream o;
        if (req>0 && k>0)
          o << "#     Spring const.     = " << k    << " kT/AA^2" << endl
            << "#     Spring eq. dist   = " << req  << " AA" << endl;
        return o.str();
      }
  };

  /*!
   * \brief Core potentials used to construct pair potentials.
   *
   * The following are low level "core" classes or snippets for handling
   * typical pair interactions. The idea is to combine these to construct
   * a full pair potential class. The combination should be done not by
   * inheritance but rather via class members.
   */
  namespace CorePotential {

    class base {
      protected:
        std::ostringstream pad;
        void setPad(char width=14, char indent=3) {
          pad << "#   " << setw(indent) << " " << setw(width);  
        }
      public:
        string name; //!< Short (preferably one-word) description of the core potential
    };

    class _hardsphere : public base {
      private:
        double inf;
      public:
        _hardsphere(double infinity=1e9) {
          name="hardsphere";
          inf=infinity;
        }
        inline double u_hs(const double &r2, double mindist) const {
          return (mindist*mindist>r2) ? 0 : inf;
        }
    };

    class lennardjones : public base {
      public:
        double eps;
        lennardjones(inputfile &in) {
          name="LJ";
          eps=in.getflt("LJeps",0.2);
        }
        inline double u_r6(const double &r2, double sigma) const {
          double x=sigma*sigma/r2;  // 2
          return x*x*x;             // 6
        }
        inline double u_r12(const double &r2, double sigma) const {
          double r6=u_r6(r2,sigma);
          return r6*r6;
        }
        inline double u_lj(const double &r2, double sigma) const {
          double r6=u_r6(r2,sigma);
          return r6*r6 - r6;
        }
    };

    class squarewell : public base {
      public:
        double threshold; //!< Threshold between particle *surface* [A]
        double depth;     //!< Energy depth [kT]
        squarewell(inputfile &in, string prefix="squarewell") {
          name="Square Well";
          threshold = in.getflt(prefix+"_threshold", 0);
          depth     = in.getflt(prefix+"_depth", 0);
        }
        inline double u_squarewell(const double &r, const double &radius1, const double &radius2) const {
          return (r-radius1-radius2<threshold) ? depth : 0;
        }
    };

    class coulomb : public base {
      public:
        double lB; //!< Bjerrum length [A]
        coulomb(inputfile &in) {
          name="Coulomb";
          lB=pyc.lB( in.getflt("epsilon_r",80.) );
        }
        inline double u_coulomb(const double &r, const double &z1, const double &z2) const { return z1*z2/r; }
    };

    class debyehuckel : public coulomb {
      protected:
        double c,k;
      public:
        debyehuckel(inputfile &in) : coulomb(in) {
          double I;
          const double zero=1e-10;
          name="Debye-Huckel";
          c=8*lB*pyc.pi*pyc.Nav/1e27;
          I=in.getflt("ionicstrength",0);  // [mol/l]
          k=sqrt( I*c );
          if (k<zero)
            k=1/in.getflt("debyelength", 1/zero); // [A]
          k=-k;
        }
        double getIonicStrength() const { return k*k/c; } //!< in [mol/l]
        double getDebyeLength() const { return -1/k; }    //!< in [A]
        inline double u_dh(const double &r, const double &z1, const double &z2) const {
          return z1*z2/r * exp(k*r);
        }
    };
  }//namespace

}//namespace
#endif
