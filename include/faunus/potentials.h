#ifndef FAU_POTENTIAL_H
#define FAU_POTENTIAL_H

#include <faunus/common.h>
#include <faunus/point.h>
#include <faunus/geometry.h>
#include <faunus/textio.h>


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

      class potbase {
        public:
          string name; //!< Short (preferably one-word) description of the core potential
          double tokT;
      };

      class hardsphere : public potbase {
        private:
          double inf;
        public:
          hardsphere(double=1e14);
          inline double u_hs(const double &r2, const double &mindist) const {
            return (mindist*mindist>r2) ? 0 : inf;
          }
      };

      class lennardjones : public potbase {
        protected:
          double eps;
        public:
          lennardjones();
          inline double r6(const double &r2, const double &sigma) const {
            double x=sigma*sigma/r2;  // 2
            return x*x*x;             // 6
          }
          inline double r12(const double &r2, const double &sigma) const {
            double x=r6(r2,sigma);
            return x*x;               // 12
          }
          inline double u_lj(const double &r2, const double &sigma, const double &eps) const {
            double x=r6(r2,sigma);
            return eps*(x*x - x);
          }
          string info(char);
      };

      class squarewell : public potbase {
        public:
          double threshold; //!< Threshold between particle *surface* [A]
          double depth;     //!< Energy depth [kT]
          squarewell(inputfile &, string="squarewell");
          inline double u_squarewell(const double &r, const double &radius1, const double &radius2) const {
            return (r-radius1-radius2<threshold) ? depth : 0;
          }

          string info(char);
      };

      class coulomb : public potbase {
        protected:
          double lB; //!< Bjerrum length [A]
        public:
          coulomb(inputfile &);
          inline double u(const double &r, const double &zz) const { return zz/r; }
          string info(char);
      };

      class debyehuckel : public coulomb {
        protected:
          double c,k;
        public:
          debyehuckel(inputfile&);
          double ionic_strength() const;
          double debye_length() const;
          inline double u(const double &r, const double &zz) const {
            return zz/r * exp(k*r);
          }
          string info(char);
      };

    } //end of Core namespace

    template<class Tgeometry, class Tcoulomb=Core::coulomb>
      class coulomb_lj {
        protected:
          Core::lennardjones lj;
          Tcoulomb el;
          const double eps;
        public:
          double tokT;
          coulomb_lj(inputfile &in) : el(in), eps(4*in.getflt("lj_eps",0.04)/el.tokT) {
            tokT=el.tokT;
          }
          inline double pairpot(const particle &a, const particle &b) const {
            double r2=Tgeometry::sqdist(a,b);
            return el.u(sqrt(r2), a.charge*b.charge) + lj.u_lj(r2, a.radius+b.radius, eps);
          }
          string info(char w=20) {
            std::ostringstream o;
            o << header("Pair Potential: " + lj.name + "+" + el.name)
              << el.info(w)
              << pad(SUB,w,"LJ epsilon") << eps*tokT << endl;
            return o.str();
          }
      };

    template<class Tgeometry, class Tcoulomb=Core::coulomb>
      class coulomb_hs {
        protected:
          Tcoulomb el;
        public:
          const double infty;
          double tokT;
          coulomb_hs(inputfile &in) : el(in), infty(1e10) { tokT=el.tokT; }
          inline double pairpot(const particle &a, const particle &b) const {
            double r2=Tgeometry::sqdist(a,b), s=a.radius+b.radius;
            if (r2<s*s)
              return infty;
            return el.u(sqrt(r2), a.charge*b.charge);
          }
          string info(char w=20) {
            std::ostringstream o;
            o << header("Pair potential: Hardsphere + " + el.name )
              << el.info(w);
            return o.str();
          }
      };

  } //end of Potential namespace

} //end of Faunus namespace
#endif
