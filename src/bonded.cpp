#include <faunus/bonded.h>

namespace Faunus {
  namespace Energy {

    HarmonicBond::HarmonicBond(double forceconst, double eqdist) : k(forceconst), req(eqdist) {}

    double HarmonicBond::energy(double r) const {
      double d=r-req;
      return k*d*d;
    }

    string HarmonicBond::info() const {
      using namespace Faunus::textio;
      std::ostringstream o;
      o << "Harmonic: k=" << k << kT << "/" << angstrom << squared << " req=" << req << _angstrom; 
      return o.str();
    }

  }//namespace
}//namespace

