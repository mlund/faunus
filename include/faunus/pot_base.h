#ifndef FAU_POT_BASE_H
#define FAU_POT_BASE_H
#include "faunus/common.h"
#include "faunus/point.h"
#include "faunus/inputfile.h"
#include "faunus/physconst.h"
#include <faunus/container.h>

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

      class base {
        public:
          void pad(std::ostringstream& o, char width) { o << "#   " << setw(width) << std::left; }
          string name; //!< Short (preferably one-word) description of the core potential
      };

      class hardsphere : public base {
        private:
          double inf;
        public:
          hardsphere(double infinity=1e9) {
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
            name="Lennard-Jones";
            eps=4*in.getflt("lj_eps",0.04);
          }
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
          string info(char w) {
            std::ostringstream o;
            pad(o,w); o << "LJ epsilon (kT)" << eps << endl;
            return o.str();
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
          string info(char w) {
            std::ostringstream o;
            pad(o,w); o << "Threshold (A)" << threshold << endl;
            pad(o,w); o << "Depth (kT)" << depth << endl;
            return o.str();
          }
       };

      class coulomb : public base {
        public:
          double lB; //!< Bjerrum length [A]
          coulomb(inputfile &in) {
            name="Coulomb";
            lB=pc::lB( in.getflt("epsilon_r",80.) );
          }
          inline double u_coulomb(const double &r, const double &zz) const { return zz/r; }

          string info(char w) {
            std::ostringstream o;
            pad(o,w); o << "Bjerrum length (A)" << lB << endl;
            return o.str();
          }
      };

      class debyehuckel : public coulomb {
        protected:
          double c,k;
        public:
          debyehuckel(inputfile &in) : coulomb(in) {
            double I;
            const double zero=1e-10;
            name="Debye-Huckel";
            c=8 * lB * pc::pi * pc::Nav / 1e27;
            I=in.getflt("dh_ionicstrength",0);  // [mol/l]
            k=sqrt( I*c );
            if (k<zero)
              k=1/in.getflt("dh_debyelength", 1/zero); // [A]
            k=-k;
          }
          double ionic_strength() const { return k*k/c; } //!< in [mol/l]
          double debye_length() const { return -1/k; }    //!< in [A]
          inline double u_dh(const double &r, const double &zz) const {
            return zz/r * exp(k*r);
          }
          string info(char w) {
            std::ostringstream o;
            o << coulomb::info(w);
            pad(o,w); o << "Ionic strength (M)" << ionic_strength() << endl;
            pad(o,w); o << "Debye length (A)" << ionic_strength() << endl;
            return o.str();
          }
      };

      class distancebase {
        protected:
          string name;
          void pad(std::ostringstream& o, char w) {
            o << "#   " << setw(w) << std::left;
          }
        public:
          virtual string info(char w) {
            std::ostringstream o;
            pad(o,w); o << "Boundaries" << name << endl;
            return o.str();
          }
      };

    } //Core namespace

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

  } //Potential namespace

} //Faunus namespace
#endif
