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

      Potbase::Potbase() { setScale(1.0); }
      void Potbase::setScale(double s)  { _setScale(s); }
      void Potbase::_setScale(double s) { tokT=s; }
      double Potbase::scale() { return tokT; }

      string Potbase::brief() {
        return _brief();
      }
      string Potbase::_brief() {
        return string("This is potential base!");
      }

      Harmonic::Harmonic(double forceconst, double eqdist) : k(forceconst), req(eqdist) {}

      void Harmonic::_setScale(double s) {
        tokT=s;
        k=k/tokT;
      }

      double Harmonic::energy(const particle &a, const particle &b, double r2) const {
        double d=sqrt(r2)-req;
        return k*d*d;
      }

      string Harmonic::_brief() {
        using namespace Faunus::textio;
        std::ostringstream o;
        o << "Harmonic: k=" << k*tokT << kT << "/" << angstrom << squared << " req=" << req << _angstrom; 
        return o.str();
      }

      HardSphere::HardSphere(double infinity) {
        name="Hardsphere";
        inf=infinity;
      }

      LennardJones::LennardJones() {
        name="Lennard-Jones";
        tokT=1;
      }

      string LennardJones::info(char w) {
        std::ostringstream o;
        return o.str();
      }

      SquareWell::SquareWell(InputMap &in, string prefix) {
        name="Square Well";
        threshold = in.get<double>(prefix+"_threshold", 0);
        depth     = in.get<double>(prefix+"_depth", 0);
      }

      string SquareWell::info(char w) {
        std::ostringstream o;
        o << pad(SUB,w,"Threshold") << threshold << " "+angstrom << endl;
        o << pad(SUB,w,"Depth") << depth << kT << endl;
        return o.str();
      }

      Coulomb::Coulomb(InputMap &in) {
        name="Coulomb";
        temp=in.get("temperature", 298.15);
        epsilon_r=in.get("epsilon_r",80.);
        pc::T = temp;
        lB=pc::lB( epsilon_r );
        tokT=lB;
      }

      string Coulomb::info(char w) {
        std::ostringstream o;
        o << pad(SUB,w,"Temperature") << temp << " K" << endl
          << pad(SUB,w,"Dielectric constant") << epsilon_r << endl
          << pad(SUB,w,"Bjerrum length") << lB << " "+angstrom << endl;
        return o.str();
      }

      DebyeHuckel::DebyeHuckel(InputMap &in) : Coulomb(in) {
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

      double DebyeHuckel::ionicStrength() const {
        return k*k/c;
      } //!< in [mol/l]

      double DebyeHuckel::debyeLength() const {
        return -1/k;
      } //!< in [A]

      string DebyeHuckel::info(char w) {
        std::ostringstream o;
        o << Coulomb::info(w);
        o << pad(SUB,w,"Ionic strength") << ionicStrength() << " mol/l" << endl;
        o << pad(SUB,w,"Debye length, 1/\u03BA") << debyeLength() << " "+angstrom << endl;
        return o.str();
      }

    } //Core namespace

  } //Potential namespace

} //Faunus namespace
#endif
