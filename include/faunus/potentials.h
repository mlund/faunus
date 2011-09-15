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
          inline double u_coulomb(const double &r, const double &zz) const { return zz/r; }
          string info(char);
      };

      class debyehuckel : public coulomb {
        protected:
          double c,k;
        public:
          debyehuckel(inputfile&);
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
        protected:
          Core::lennardjones lj;
          Core::coulomb el;
          const double eps;
        public:
          double tokT;
          coulomb_lj(inputfile &in) : el(in), eps(4*in.getflt("lj_eps",0.04)/el.tokT) {
            tokT=el.tokT;
          }
          inline double pairpot(const particle &a, const particle &b) const {
            double r2=Tgeometry::sqdist(a,b);
            return el.u_coulomb(sqrt(r2), a.charge*b.charge)
              + lj.u_lj(r2, a.radius+b.radius, eps);
          }
          string info(char w=20) {
            std::ostringstream o;
            o << "\n# PAIR POTENTIAL: " << lj.name << "+" << el.name << "\n" << el.info(w);
            el.pad(o,w); o << "LJ epsilon" << eps*tokT << endl;
            return o.str();
          }
      };

  } //end of Potential namespace

} //end of Faunus namespace
#endif
