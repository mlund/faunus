#include <faunus/common.h>
#include <faunus/point.h>
#include <faunus/group.h>
#include <faunus/energy.h>
#include <faunus/textio.h>

namespace Faunus {
  namespace Energy {
    Energybase::Energybase() : geo(NULL) {}

    Energybase::~Energybase() {}

    Geometry::Geometrybase& Energybase::getGeometry() {
      return *geo;
    }

    // external particle interactions
    double Energybase::all2p(const p_vec &p, const particle &a) {return 0;}
    double Energybase::p2p(const particle &a, const particle &b) {return 0;}

    // single particle interactions
    double Energybase::all2all(const p_vec &p) {return 0;}
    double Energybase::i2i(const p_vec &p, int i, int j) {return 0;}
    double Energybase::i2g(const p_vec &p, const group &g, int i) {return 0;}
    double Energybase::i2all(const p_vec &p, int i) {return 0;}
    double Energybase::i_external(const p_vec &p, int i) {return 0;}
    double Energybase::i_internal(const p_vec &p, int i) {return 0;}

    // group interactions
    double Energybase::g2g(const p_vec &p, const group &g1, const group &g2) {return 0;}
    double Energybase::g2all(const p_vec &p, const group &g) {return 0;}
    double Energybase::g_external(const p_vec &p, const group &g) {return 0;}
    double Energybase::g_internal(const p_vec &p, const group &g) {return 0;}
    string Energybase::info() {
      std::ostringstream o;
      return o.str();
    }

  }//namespace
}//namespace

//dynamic cast of group->derive class
//const molecular* m = dynamic_cast<const molecular*>(&g1);


