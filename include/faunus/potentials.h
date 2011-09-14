#ifndef FAU_POTENTIAL_H
#define FAU_POTENTIAL_H

#include <faunus/common.h>
#include <faunus/point.h>
#include <faunus/geometry.h>

namespace Faunus {

  namespace Potential {

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
          void pad(std::ostringstream&, char);
          string name; //!< Short (preferably one-word) description of the core potential
          double tokT;
      };

      class hardsphere : public potbase {
        private:
          double inf;
        public:
          hardsphere(double=1e9);
          inline double u_hs(const double &r2, double mindist) const {
            return (mindist*mindist>r2) ? 0 : inf;
          }
      };

      class lennardjones : public potbase {
        public:
          double eps;
          lennardjones(inputfile &);
          inline double u_r6(const double &r2, const double &sigma) const {
            double x=sigma*sigma/r2;  // 2
            return x*x*x;             // 6
          }

          inline double u_r12(const double &r2, const double &sigma) const {
            double r6=u_r6(r2,sigma);
            return r6*r6;
          }

          inline double u_lj(const double &r2, const double &sigma) const {
            double r6=u_r6(r2,sigma);
            return r6*r6 - r6;
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
        public:
          double lB; //!< Bjerrum length [A]
          coulomb(inputfile &);
          inline double u_coulomb(const double &r, const double &zz) const { return zz/r; }
          string info(char);
      };

      class debyehuckel : public coulomb {
        protected:
          double c,k;
        public:
          debyehuckel(inputfile &);
          double ionic_strength() const;
          double debye_length() const;
          inline double u_dh(const double &r, const double &zz) const {
            return zz/r * exp(k*r);
          }
          string info(char);
      };

    } //end of Core namespace

    template<class Tgeometry>
      class coulomb_lj {
        public:
          Core::lennardjones lj;
          Core::coulomb el;
          coulomb_lj(inputfile &in) : lj(in), el(in) {}
          inline double pairpot(const particle &a, const particle &b) const {
            double r2=a.sqdist<Tgeometry>(b);
            return el.u_coulomb(sqrt(r2), a.charge*b.charge) + lj.u_lj(r2,a.radius+b.radius);
          }
          string info(char w=20) {
            return "\n# PAIR POTENTIAL: " + lj.name +"+" + el.name + "\n"
              + lj.info(w) + el.info(w);
          }
      };

  } //end of Potential namespace

} //end of Faunus namespace
#endif
