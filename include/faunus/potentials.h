#ifndef FAU_POTENTIAL_H
#define FAU_POTENTIAL_H

#include <faunus/common.h>
#include <faunus/point.h>
#include <faunus/geometry.h>
#include <faunus/textio.h>
#include <faunus/inputfile.h>

namespace Faunus {

  namespace Potential {

    using namespace Faunus::textio;

    /*!
     * \brief Core potentials used to construct pair potentials.
     *
     * The following are low level "core" classes or snippets for handling
     * typical pair interactions. The idea is to combine these to construct
     * a full pair potential class. The combination should be done not by
     * inheritance but rather via class members.
     */
    namespace Core {

      class Potbase {
        private:
          virtual string _brief();
          virtual void _setScale(double);
        public:
          double tokT; //!< to be removed!
          string name; //!< Short (preferably one-word) description of the core potential
          Potbase();
          string brief();
          virtual double energy(const particle&, const particle&, double) const=0;
          void setScale(double=1);
          double scale();
      };

      class Harmonic : public Potbase {
        private:
          string _brief();
          void _setScale(double);
          double _scale();
        public:
          double k;   //!< Force constant (kT/AA^2)
          double req; //!< Equilibrium distance (AA)
          Harmonic(double=0, double=0);
          double energy(const particle&, const particle&, double) const;
      };

      class HardSphere : public Potbase {
        private:
          double inf;
        public:
          HardSphere(double=1e14);
          inline double energy(const particle &a, const particle &b, double r2) const {
            double mindist=a.radius+b.radius;
            if (mindist*mindist<r2)
              return inf;
            return 0;
          }
      };

      class LennardJones : public Potbase {
        protected:
          double eps;
        public:
          LennardJones();
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
          inline double energy(const particle &a, const particle &b, double r2) const {
            return energy(a.radius+b.radius, r2);
          }
          string info(char);
      };

      class SquareWell : public Potbase {
        public:
          double threshold; //!< Threshold between particle *surface* [A]
          double depth;     //!< Energy depth [kT]
          SquareWell(InputMap&, string="SquareWell");
          inline double energy(const particle &a, const particle &b, double r2) const {
            return ( sqrt(r2)-a.radius-b.radius<threshold ) ? depth : 0;
          }

          string info(char);
      };

      class Coulomb : public Potbase {
        private:
          double temp, epsilon_r;
        protected:
          double lB; //!< Bjerrum length [A]
        public:
          Coulomb(InputMap &);
          inline double energy(double zz, double r) const { return zz/r; } 
          inline double energy(const particle &a, const particle &b, double r2) const {
            return energy( a.charge*b.charge, sqrt(r2) );
          }
          string info(char);
      };

      class DebyeHuckel : public Coulomb {
        protected:
          double c,k;
        public:
          DebyeHuckel(InputMap&);
          double ionicStrength() const;
          double debyeLength() const;
          inline double energy(double zz, double r) const { return zz/r * exp(k*r); }
          inline double energy(const particle &a, const particle &b, double r2) const {
            r2=sqrt(r2);
            return a.charge*b.charge/r2 * exp(k*r2);
          }
          string info(char);
      };

    } //end of Core namespace

    template<class Tgeometry, class Tcoulomb=Core::Coulomb>
      class CoulombLJ {
        protected:
          Core::LennardJones lj;
          Tcoulomb el;
          const double eps;
        public:
          string name; //!< Single line describing the potential
          Tgeometry geo;
          double tokT;
          CoulombLJ(InputMap &in) : el(in), eps(4*in.get("lj_eps",0.04)/el.tokT), geo(in) { 
            tokT=el.tokT;
            name=lj.name+"+"+el.name;
          }
          inline double energy(const particle &a, const particle &b) const {
            double r2=geo.sqdist(a,b);
            return el.energy(a.charge*b.charge,sqrt(r2)) + lj.energy(a.radius+b.radius, r2);
          }
          inline double excluded(const particle &a, const particle &b) const {
            return lj.energy( a.radius+b.radius, geo.sqdist(a,b) );
          }
          string info(char w=20) {
            std::ostringstream o;
            o << pad(SUB,w,"Pair potential:") << name << endl
              << el.info(w)
              << pad(SUB,w,"LJ epsilon") << eps*tokT << kT << endl;
            return o.str();
          }
      };

    template<class Tgeometry, class Tcoulomb=Core::Coulomb>
      class CoulombHS {
        protected:
          Tcoulomb el;
        public:
          string name; //!< Single line describing the potential
          const double infty;
          double tokT;
          CoulombHS(InputMap &in) : el(in), infty(1e10) {
            name="Hardsphere + " + el.name;
            tokT=el.tokT;
          }
          inline double energy(const particle &a, const particle &b) {
            double r2=Tgeometry::sqdist(a,b), s=a.radius+b.radius;
            if (r2<s*s)
              return infty;
            return el.energy( a.charge*b.charge, sqrt(r2) );
          }
          inline double excluded(const particle &a, const particle &b) {
            double r2=Tgeometry::sqdist(a,b), s=a.radius+b.radius;
            if (r2<s*s)
              return infty;
            return 0;
          }
          string info(char w=20) {
            std::ostringstream o;
            o << pad(SUB,w,"Pair potential:") << name << endl
              << el.info(w);
            return o.str();
          }
      };

  } //end of Potential namespace

} //end of Faunus namespace
#endif
