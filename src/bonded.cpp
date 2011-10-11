#include <faunus/bonded.h>
#include <faunus/geometry.h>
#include <faunus/textio.h>

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

    ParticleBonds::ParticleBonds() {
      pairs::name+="Particle Bonds";
    }

    double ParticleBonds::bondEnergy(Geometry::Geometrybase &geo, const p_vec &p, int i) {
      double u=0;
      for (auto &j : pairs::list[i])
        u+=j.second->energy( geo.dist( p[i], p[j.first] ) );
      return u;
    }

  }//namespace
}//namespace

