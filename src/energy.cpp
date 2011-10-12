#include <faunus/common.h>
#include <faunus/point.h>
#include <faunus/group.h>
#include <faunus/energy.h>
#include <faunus/textio.h>

namespace Faunus {
  namespace Energy {
    Energybase::Energybase() : geo(NULL) {
      name="Hamiltonian";
    }

    Energybase::~Energybase() {}

    void Energybase::setGeometry(Geometry::Geometrybase &g) {
      geo=&g;
    }

    Geometry::Geometrybase& Energybase::getGeometry() {
      return *geo;
    }

    // external particle interactions
    double Energybase::all2p(const p_vec &p, const particle &a) { return 0; }
    double Energybase::p2p(const particle &a, const particle &b) { return 0; }
    double Energybase::v2v(const p_vec &v1, const p_vec &v2) { return 0; }
    double Energybase::p_external(const particle &a) { return 0; }

    // single particle interactions
    double Energybase::all2all(const p_vec &p) {return 0;}
    double Energybase::i2i(const p_vec &p, int i, int j) {return 0;}
    double Energybase::i2g(const p_vec &p, Group &g, int i) {return 0;}
    double Energybase::i2all(const p_vec &p, int i) {return 0;}
    double Energybase::i_external(const p_vec &p, int i) {return 0;}
    double Energybase::i_internal(const p_vec &p, int i) {return 0;}
    double Energybase::i_total(const p_vec &p, int i) {
      return i2all(p,i) + i_external(p,i) + i_internal(p,i);
    }

    // Group interactions
    double Energybase::g2g(const p_vec &p, Group &g1, Group &g2) {return 0;}
    double Energybase::g2all(const p_vec &p, Group &g) {return 0;}
    double Energybase::g_external(const p_vec &p, Group &g) {return 0;}
    double Energybase::g_internal(const p_vec &p, Group &g) {return 0;}

    string Energybase::info() {
      using namespace textio;
      std::ostringstream o;
      o << header("Energy: " + name); 
      return o.str();
    }

    Bonded::Bonded() {
      name="Bonded";
      geo=NULL;
    }

    Bonded::Bonded(Geometry::Geometrybase &g) {
      name="Bonded";
      geo=&g;
    }

    double Bonded::i_internal(const p_vec &p, int i) {
      return bonds.totalEnergy(*geo, p, i);
    }

    double Bonded::g_internal(const p_vec &p, Group &g) {
      return bonds.totalEnergy(*geo, p, g);
    }


    /*!
     * This adds an Energybase class to the Hamiltonian. If the geometry of the
     * added class is undefined it will get the current geometry of the Hamiltonian.
     * If the geometry of the added energybase is defined, this geometry is copied to
     * the Hamiltonian.
     */
    void Hamiltonian::add(Energybase &e) {
      assert(&e!=NULL);
      if (&e.getGeometry()!=NULL)
        geo=&e.getGeometry();
      else
        e.setGeometry(*geo);
      baselist.push_back(&e);
    }

    double Hamiltonian::all2all(const p_vec &p) {
      double u=0;
      for (auto &b : baselist)
        u += b->all2all(p);
      return u;
    }

    double Hamiltonian::p2p(const particle &p1, const particle &p2) {
      double u=0;
      for (auto &b : baselist)
        u += b->p2p( p1,p2 );
      return u;
    }

    double Hamiltonian::all2p(const p_vec &p, const particle &a) {
      double u=0;
      for (auto &b : baselist)
        u += b->all2p( p,a );
      return u;
    }

    // single particle interactions
    double Hamiltonian::i2i(const p_vec &p, int i, int j) {
      double u=0;
      for (auto &b : baselist)
        u += b->i2i( p,i,j );
      return u;
    }

    double Hamiltonian::i2g(const p_vec &p, Group &g, int i) {
      double u=0;
      for (auto &b : baselist)
        u += b->i2g(p,g,i);
      return u;
    }

    double Hamiltonian::i2all(const p_vec &p, int i) {
      double u=0;
      for (auto &b : baselist)
        u += b->i2all(p,i);
      return u;
    }

    double Hamiltonian::i_external(const p_vec &p, int i) {
      double u=0;
      for (auto &b : baselist)
        u += b->i_external(p,i);
      return u;
    }
    double Hamiltonian::i_internal(const p_vec &p, int i) {
      double u=0;
      for (auto &b : baselist)
        u += b->i_internal(p,i);
      return u;
    }

    // Group interactions
    double Hamiltonian::g2g(const p_vec &p, Group &g1, Group &g2) {
      double u=0;
      for (auto &b : baselist)
        u += b->g2g(p,g1,g2);
      return u;
    }

    double Hamiltonian::g2all(const p_vec &p, Group &g) {
      double u=0;
      for (auto &b : baselist)
        u += b->g2all(p,g);
      return u;
    }

    double Hamiltonian::g_external(const p_vec &p, Group &g) {
      double u=0;
      for (auto &b : baselist)
        u += b->g_external(p,g);
      return u;
    }

    double Hamiltonian::g_internal(const p_vec &p, Group &g) {
      double u=0;
      for (auto &b : baselist)
        u += b->g_internal(p,g);
      return u;
    }

    string Hamiltonian::info() {
      using namespace textio;
      char w=25;
      std::ostringstream o;
      o << Energybase::info();
      o << indent(SUB) << "Registered Energy Functions:"<< endl;
      int i=1;
      for (auto e : baselist)
        o << indent(SUBSUB) << std::left << setw(4) << i++ << e->name << endl;
      for (auto e : baselist)
        o << e->info();
      return o.str();
    }


  }//namespace
}//namespace

//dynamic cast of Group->derive class
//const molecular* m = dynamic_cast<const molecular*>(&g1);


