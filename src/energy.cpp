#include <faunus/common.h>
#include <faunus/point.h>
#include <faunus/group.h>
#include <faunus/energy.h>

namespace Faunus {
  namespace Energy {

    Geometry::geometrybase* energybase::geo;

    // single particle interactions
    double energybase::all2all(const vector<particle> &p) {return 0;}
    double energybase::i2i(const vector<particle> &p, int i, int j) {return 0;}
    double energybase::i2g(const vector<particle> &p, const group &g, int i) {return 0;}
    double energybase::i2all(const vector<particle> &p, int i) {return 0;}
    double energybase::i_external(const vector<particle> &p, int i) {return 0;}
    double energybase::i_internal(const vector<particle> &p, int i) {return 0;}

    // group interactions
    double energybase::g2g(const vector<particle> &p, const group &g1, const group &g2) {return 0;}
    double energybase::g2all(const vector<particle> &p, const group &g) {return 0;}
    double energybase::g_external(const vector<particle> &p, const group &g) {return 0;}
    double energybase::g_internal(const vector<particle> &p, const group &g) {return 0;}
    string energybase::info() {
      std::ostringstream o;
      return o.str();
    }

  }//namespace
}//namespace

//dynamic cast of group->derive class
//const molecular* m = dynamic_cast<const molecular*>(&g1);


