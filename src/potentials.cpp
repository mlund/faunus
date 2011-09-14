#ifndef FAU_POT_BASE_H
#define FAU_POT_BASE_H

#include <faunus/common.h>
#include <faunus/point.h>
#include <faunus/inputfile.h>
#include <faunus/physconst.h>
#include <faunus/geometry.h>
#include <faunus/potentials.h>

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

      void potbase::pad(std::ostringstream& o, char width) {
        o << "#   " << setw(width) << std::left;
      }

      hardsphere::hardsphere(double infinity) {
        name="hardsphere";
        inf=infinity;
      }

      lennardjones::lennardjones(inputfile &in) {
        name="Lennard-Jones";
        eps=4*in.getflt("lj_eps",0.04);
      }

      string lennardjones::info(char w) {
        std::ostringstream o;
        pad(o,w); o << "LJ epsilon (kT)" << eps << endl;
        return o.str();
      }

      squarewell::squarewell(inputfile &in, string prefix) {
        name="Square Well";
        threshold = in.getflt(prefix+"_threshold", 0);
        depth     = in.getflt(prefix+"_depth", 0);
      }

      string squarewell::info(char w) {
        std::ostringstream o;
        pad(o,w); o << "Threshold (A)" << threshold << endl;
        pad(o,w); o << "Depth (kT)" << depth << endl;
        return o.str();
      }

      coulomb::coulomb(inputfile &in) {
        name="Coulomb";
        lB=pc::lB( in.getflt("epsilon_r",80.) );
        tokT=lB;
      }

      string coulomb::info(char w) {
        std::ostringstream o;
        pad(o,w); o << "Bjerrum length (A)" << lB << endl;
        return o.str();
      }

      debyehuckel::debyehuckel(inputfile &in) : coulomb(in) {
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

      double debyehuckel::ionic_strength() const {
        return k*k/c;
      } //!< in [mol/l]

      double debyehuckel::debye_length() const {
        return -1/k;
      } //!< in [A]

      string debyehuckel::info(char w) {
        std::ostringstream o;
        o << coulomb::info(w);
        pad(o,w); o << "Ionic strength (M)" << ionic_strength() << endl;
        pad(o,w); o << "Debye length (A)" << ionic_strength() << endl;
        return o.str();
      }

    } //Core namespace

  } //Potential namespace

} //Faunus namespace
#endif
