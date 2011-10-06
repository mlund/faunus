#ifndef FAU_POT_BASE_H
#define FAU_POT_BASE_H

#include <faunus/common.h>
#include <faunus/point.h>
#include <faunus/inputfile.h>
#include <faunus/physconst.h>
#include <faunus/geometry.h>
#include <faunus/potentials.h>
#include <faunus/textio.h>

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

      hardsphere::hardsphere(double infinity) {
        name="hardsphere";
        inf=infinity;
      }

      lennardjones::lennardjones() {
        name="Lennard-Jones";
        tokT=1;
      }

      string lennardjones::info(char w) {
        std::ostringstream o;
        return o.str();
      }

      squarewell::squarewell(InputMap &in, string prefix) {
        name="Square Well";
        threshold = in.get<double>(prefix+"_threshold", 0);
        depth     = in.get<double>(prefix+"_depth", 0);
      }

      string squarewell::info(char w) {
        std::ostringstream o;
        o << pad(SUB,w,"Threshold") << threshold << " "+angstrom << endl;
        o << pad(SUB,w,"Depth") << depth << kT << endl;
        return o.str();
      }

      coulomb::coulomb(InputMap &in) {
        name="Coulomb";
        pc::T=in.get<double>("temperature", 298.15);
        lB=pc::lB( in.get<double>("epsilon_r",80.) );
        tokT=lB;
      }

      string coulomb::info(char w) {
        std::ostringstream o;
        o << pad(SUB,w,"Bjerrum length") << lB << " "+angstrom << endl;
        return o.str();
      }

      debyehuckel::debyehuckel(InputMap &in) : coulomb(in) {
        double I;
        const double zero=1e-10;
        name="Debye-Huckel";
        c=8 * lB * pc::pi * pc::Nav / 1e27;
        I=in.get<double>("dh_ionicstrength",0);  // [mol/l]
        k=sqrt( I*c );
        if (k<zero)
          k=1/in.get<double>("dh_debyelength", 1/zero); // [A]
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
        o << pad(SUB,w,"Ionic strength") << ionic_strength() << " mol/l" << endl;
        o << pad(SUB,w,"Debye length, 1/\u03BA") << debye_length() << " "+angstrom << endl;
        return o.str();
      }

    } //Core namespace

  } //Potential namespace

} //Faunus namespace
#endif
