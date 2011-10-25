#ifndef FAU_POTENTIAL_H
#define FAU_POTENTIAL_H

#include <faunus/common.h>
#include <faunus/point.h>
#include <faunus/geometry.h>
#include <faunus/textio.h>
#include <faunus/inputfile.h>
#include <faunus/physconst.h>

namespace Faunus {

  namespace Potential {

    class DebyeHuckel;

    using namespace Faunus::textio;

    class PairPotentialBase {
      private:
        virtual string _brief()=0;
        virtual void _setScale(double);
      protected:
        double _tokT;
      public:  
        string name; //!< Short (preferably one-word) description of the core potential
        PairPotentialBase();
        string brief();
        void setScale(double=1);
        double tokT();
        virtual double operator() (const particle&, const particle&, double) const=0;
    };

    class Harmonic : public PairPotentialBase {
      private:
        string _brief();
        void _setScale(double);
      public:
        double k;   //!< Force constant (kT/AA^2)
        double req; //!< Equilibrium distance (AA)
        Harmonic(double=0, double=0);
        double operator() (const particle&, const particle&, double) const; //!< Pair interaction energy (kT)
    };

    class HardSphere : public PairPotentialBase {
      private:
        string _brief();
      public:
        HardSphere();
        HardSphere(InputMap&);
        inline double operator() (const particle &a, const particle &b, double r2) const {
          double mindist=a.radius+b.radius;
          if (r2<mindist*mindist)
            return pc::infty;
          return 0;
        }
        string info(char w);
    };

    class LennardJones : public PairPotentialBase {
      private:
        string _brief();
        void _setScale(double);
      protected:
        double eps;
      public:
        LennardJones();
        LennardJones(InputMap&);
        inline double r6(double &sigma, double r2) const {
          double x=sigma*sigma/r2;  // 2
          return x*x*x;             // 6
        }
        inline double r12(double sigma, double r2) const {
          double x=r6(sigma,r2);
          return x*x;               // 12
        }
        inline double energy(double sigma, double r2) const {
          double x=r6(sigma,r2);
          return eps*(x*x - x);
        }
        inline double operator() (const particle &a, const particle &b, double r2) const {
          return energy(a.radius+b.radius, r2);
        }
        string info(char);
    };

    class SquareWell : public PairPotentialBase {
      private:
        string _brief();
        void _setScale(double);
      public:
        double threshold; //!< Threshold between particle *surface* [A]
        double depth;     //!< Energy depth [kT]
        SquareWell(InputMap&, string="SquareWell");
        inline double operator() (const particle &a, const particle &b, double r2) const {
          return ( sqrt(r2)-a.radius-b.radius<threshold ) ? depth : 0;
        }
        string info(char);
    };

    class Coulomb : public PairPotentialBase {
      friend class DebyeHuckel;
      private:
      string _brief();
      void _setScale(double);
      double temp, epsilon_r;
      protected:
      double lB; //!< Bjerrum length [A]
      public:
      Coulomb(InputMap &);
      inline double energy(double zz, double r) const { return zz/r; } 
      inline double operator() (const particle &a, const particle &b, double r2) const {
        return energy( a.charge*b.charge, sqrt(r2) );
      }
      string info(char);
    };

    class DebyeHuckel : public Coulomb {
      private:
        string _brief();
      protected:
        double c,k;
      public:
        DebyeHuckel(InputMap&);
        double ionicStrength() const;
        double debyeLength() const;
        inline double energy(double zz, double r) const { return zz/r * exp(k*r); }
        inline double operator() (const particle &a, const particle &b, double r2) const {
          r2=sqrt(r2);
          return a.charge*b.charge/r2 * exp(k*r2);
        }
        string info(char);
    };


    template<class Tgeometry, class Tcoulomb=Coulomb, class Tshortranged=LennardJones>
      class CoulombSR : public PairPotentialBase {
        private:
          string _brief() { return name; };
        protected:
          Tshortranged sr;
          Tcoulomb el;
        public:
          Tgeometry geo;
          CoulombSR(InputMap &in) : sr(in), el(in), geo(in) {
            setScale( el.tokT() );
            sr.setScale( el.tokT() );
            name=sr.name+"+"+el.name;
          }
          inline double energy(const particle &a, const particle &b) const {
            double r2=geo.sqdist(a,b);
            return el.energy( a.charge*b.charge, sqrt(r2) ) + sr(a,b,r2);
          }
          inline double operator() (const particle &a, const particle &b, double r2) const {
            return el.energy( a.charge*b.charge, sqrt(r2) ) + sr(a,b,r2);
          }
          string info(char w=20) {
            std::ostringstream o;
	    o << el.info(w) << sr.info(w);
            return o.str();
          }
      };
      
  } //end of potential namespace

} //end of Faunus namespace
#endif
