#include <faunus/bonded.h>
#include <faunus/geometry.h>
#include <faunus/textio.h>
#include <faunus/group.h>
#include <faunus/potentials.h>

namespace Faunus {
  namespace Energy {

    Bondbase::~Bondbase() { }

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

    double ParticleBonds::totalEnergy(Geometry::Geometrybase &geo, const p_vec &p, int i) {
      assert( &geo!=NULL );  //debug
      assert( i<(int)p.size() ); //debug
      double u=0;
      for (auto &m2 : pairs::list[i])
        u+=m2.second->energy( p[i], p[m2.first], geo.sqdist( p[i], p[m2.first] ) );
      return u;
    }

    //!< Return the total bond energy by considering each bond pair exactly once (kT)
    double ParticleBonds::totalEnergy(Geometry::Geometrybase &geo, const p_vec &p) {
      assert(&geo!=NULL);  //debug
      std::set<int> done;
      double u=0;
      for (auto &m1 : pairs::list) {
        for (auto &m2 : m1.second ) {
          if ( done.find(m2.first)==done.end() ) {
            assert(m1.first<(int)p.size()); //debug
            assert(m2.first<(int)p.size()); //debug
            u += m2.second->energy( p[m1.first], p[m2.first], geo.sqdist( p[m1.first], p[m2.first] ) );
          }
        }
        done.insert(m1.first); // exclude index from subsequent loops
      }
      return u;
    }

    double ParticleBonds::totalEnergy(Geometry::Geometrybase &geo, const p_vec &p, const Group &g) {
      assert(&geo!=NULL);  //debug
      std::set<int> done;
      double u=0;
      for (auto i=g.beg; i<=g.end; i++) {
        for (auto &m2 : pairs::list[i]) {
          if ( done.find(m2.first)==done.end() ) {
            u += m2.second->energy( p[i], p[m2.first], geo.sqdist( p[i], p[m2.first] ) );
          }
        }
        done.insert(i);
      }
      return u;
    }

  }//namespace
}//namespace

