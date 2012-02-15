#include <faunus/common.h>
#include <faunus/point.h>
#include <faunus/energy.h>
#include <faunus/textio.h>

namespace Faunus {
  namespace Energy {
    Energybase::Energybase() : geo(nullptr) {
      w=25;
    }

    Energybase::~Energybase() {
    }

    bool Energybase::setGeometry(Geometry::Geometrybase &g) {
      geo=&g;
      return true;
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

    double Energybase::external() {return 0;}

    string Energybase::info() {
      assert(!name.empty() && "Assign a name to energy class!");
      std::ostringstream o;
      if (!_info().empty())
        o << textio::header("Energy: " + name)
          << _info();
      return o.str();
    }

    Bonded::Bonded() {
      name="Bonded Particles";
      geo=nullptr;
    }

    Bonded::Bonded(Geometry::Geometrybase &g) {
      name="Bonded";
      geo=&g;
    }

    double Bonded::i2all(const p_vec &p, int i) {
      return bonds.totalEnergy(*geo, p, i);
    }

    double Bonded::g_internal(const p_vec &p, Group &g) {
      return bonds.totalEnergy(*geo, p, g);
    }

    string Bonded::_info() {
      return bonds.info();
    }

    ExternalPressure::ExternalPressure(Geometry::Geometrybase &e, double pressure) {
      assert(&e!=nullptr && "Geometry must be defined for this energy term");
      assert(pressure>=0 && "Specify non-negative pressure");
      name="External Pressure";
      setGeometry(e);
      P=pressure;
    }

    double ExternalPressure::external() {
      double V=geo->getVolume();
      assert(V!=0 && P!=0 && "Pressure and/or volume is zero");
      return P*V - log(V); 
    }

    double ExternalPressure::g_external(const p_vec &p, Group &g) {
      double N=1,
             V=geo->getVolume();
      if (g.id==Group::ATOMIC)
        N=g.size();
      return -N*log(V);
    }

    string ExternalPressure::_info() {
      char w=15;
      using namespace textio;
      std::ostringstream o;
      o << pad(SUB,w,"Pressure") << P*1e30/pc::Nav << " mM = " << P*pc::kB*pc::T*1e30 << " Pa" << endl;
      return o.str();
    }

    Hamiltonian::Hamiltonian() {
      name="Hamiltonian";
    }

    /*!
     * This adds an existing Energybase child to the Hamiltonian. If the geometry of the
     * added class is undefined it will get the current geometry of the Hamiltonian.
     * If the geometry of the added energybase is defined, and the Hamiltonian is empty, the added
     * geometry is copied to the Hamiltonian.
     */
    void Hamiltonian::add(Energybase &e) {
      if (&e.getGeometry()==nullptr) {
        if (baselist.empty()) {
          string s="Error! First Energybase class must have a well defined geometry!\n";
          assert(!"Error!");
          cerr << s;
          return;
        }
        e.setGeometry(*geo);
      }
      else if (geo==nullptr)
        geo = &e.getGeometry();
      baselist.push_back( &e );
    }

    void Hamiltonian::setVolume(double V) {
      for (auto b : baselist )
        if ( &b->getGeometry() != nullptr )
          b->getGeometry().setVolume(V);
    }

    double Hamiltonian::all2all(const p_vec &p) {
      double u=0;
      for (auto b : baselist)
        u += b->all2all(p);
      return u;
    }

    double Hamiltonian::p2p(const particle &p1, const particle &p2) {
      double u=0;
      for (auto b : baselist)
        u += b->p2p( p1,p2 );
      assert(u!=0);
      return u;
    }

    double Hamiltonian::all2p(const p_vec &p, const particle &a) {
      double u=0;
      for (auto b : baselist)
        u += b->all2p( p,a );
      return u;
    }

    // single particle interactions
    double Hamiltonian::i2i(const p_vec &p, int i, int j) {
      double u=0;
      for (auto b : baselist)
        u += b->i2i( p,i,j );
      return u;
    }

    double Hamiltonian::i2g(const p_vec &p, Group &g, int i) {
      double u=0;
      for (auto b : baselist)
        u += b->i2g(p,g,i);
      return u;
    }

    double Hamiltonian::i2all(const p_vec &p, int i) {
      double u=0;
      for (auto b : baselist)
        u += b->i2all(p,i);
      return u;
    }

    double Hamiltonian::i_external(const p_vec &p, int i) {
      double u=0;
      for (auto b : baselist)
        u += b->i_external(p,i);
      return u;
    }
    double Hamiltonian::i_internal(const p_vec &p, int i) {
      double u=0;
      for (auto b : baselist)
        u += b->i_internal(p,i);
      return u;
    }

    // Group interactions
    double Hamiltonian::g2g(const p_vec &p, Group &g1, Group &g2) {
      double u=0;
      for (auto b : baselist)
        u += b->g2g(p,g1,g2);
      return u;
    }

    double Hamiltonian::g2all(const p_vec &p, Group &g) {
      double u=0;
      for (auto b : baselist)
        u += b->g2all(p,g);
      return u;
    }

    double Hamiltonian::g_external(const p_vec &p, Group &g) {
      double u=0;
      for (auto b : baselist)
        u += b->g_external(p,g);
      return u;
    }

    double Hamiltonian::g_internal(const p_vec &p, Group &g) {
      double u=0;
      for (auto b : baselist)
        u += b->g_internal(p,g);
      return u;
    }

    double Hamiltonian::external() {
      double u=0;
      for (auto b : baselist)
        u += b->external();
      return u;
    }

    double Hamiltonian::v2v(const p_vec &v1, const p_vec &v2) {
      double u=0;
      for (auto b : baselist)
        u += b->v2v(v1,v2);
      return u;
    }

    string Hamiltonian::_info() {
      using namespace textio;
      std::ostringstream o;
      o << indent(SUB) << "Registered Energy Functions:"<< endl;
      int i=1;
      for (auto e : baselist)
        o << indent(SUBSUB) << std::left << setw(4) << i++ << e->name << endl;
      for (auto e : baselist)
        o << e->info();
      return o.str();
    }

    ParticleBonds::ParticleBonds() {
      pairs::name+="Particle Bonds";
    }

    double ParticleBonds::i2i(Geometry::Geometrybase &geo, const p_vec &p, int i, int j) {
      assert(!"unimplemented!");
      return 0;
    }

    double ParticleBonds::totalEnergy(Geometry::Geometrybase &geo, const p_vec &p, int i) {
      assert( &geo!=nullptr );   //debug
      assert( i<(int)p.size() ); //debug
      double u=0;
      for (auto &m2 : pairs::list[i])
        u+=m2.second->tokT() * m2.second->operator()( p[i], p[m2.first], geo.sqdist( p[i], p[m2.first] ) );
      return u;
    }

    //!< Return the total bond energy by considering each bond pair exactly once (kT)
    double ParticleBonds::totalEnergy(Geometry::Geometrybase &geo, const p_vec &p) {
      assert(&geo!=nullptr);  //debug
      std::set<int> done;
      double u=0;
      for (auto &m1 : pairs::list) {
        for (auto &m2 : m1.second ) {
          if ( done.find(m2.first)==done.end() ) {
            assert(m1.first<(int)p.size()); //debug
            assert(m2.first<(int)p.size()); //debug
            u += m2.second->tokT() * m2.second->operator()( p[m1.first], p[m2.first], geo.sqdist( p[m1.first], p[m2.first] ) );
          }
        }
        done.insert(m1.first); // exclude index from subsequent loops
      }
      return u;
    }

    double ParticleBonds::totalEnergy(Geometry::Geometrybase &geo, const p_vec &p, const Group &g) {
      assert(&geo!=nullptr);  //debug
      std::set<int> done;
      double u=0;
      for (auto i=g.front(); i<=g.back(); i++) {
        for (auto &m2 : pairs::list[i]) {
          if ( done.find(m2.first)==done.end() ) {
            u += m2.second->tokT() * m2.second->operator()( p[i], p[m2.first], geo.sqdist( p[i], p[m2.first] ) );
          }
        }
        done.insert(i);
      }
      return u;
    }

    EnergyRest::EnergyRest() {
      usum=0;
      name="Dummy energy (drift calc. aid)";
    }

    void EnergyRest::add(double du) { usum+=du; }

    double EnergyRest::external() {
      return usum;
    }

    string EnergyRest::_info(){ return ""; }

    /*!
     * The InputMap is searched for the following keywords that defines the two points
     * that spans the allowed region:
     *
     * \li \c constrain_upper.x
     * \li \c constrain_upper.y
     * \li \c constrain_upper.z
     * \li \c constrain_lower.x
     * \li \c constrain_lower.y
     * \li \c constrain_lower.z
     *
     * The values in the upper point must be higher than those in the lower point. The default
     * values are +infinity for the upper point and -infinity for the lower, meaning that all
     * space is available.
     */
    RestrictedVolume::RestrictedVolume(InputMap &in) {
      name="Restricted Volume";
      lower.x = in.get<double>("constrain_lower.x", -pc::infty);
      lower.y = in.get<double>("constrain_lower.y", -pc::infty);
      lower.z = in.get<double>("constrain_lower.z", -pc::infty);
      upper.x = in.get<double>("constrain_upper.x", pc::infty);
      upper.y = in.get<double>("constrain_upper.y", pc::infty);
      upper.z = in.get<double>("constrain_upper.z", pc::infty);
      assert(upper.x>lower.x && "Upper bound must be bigger than lower bound!");
      assert(upper.y>lower.y && "Upper bound must be bigger than lower bound!");
      assert(upper.z>lower.z && "Upper bound must be bigger than lower bound!");
    }

    string RestrictedVolume::_info() {
      using namespace textio;
      char w=15;
      string d=" x ";
      std::ostringstream o;
      o << indent(SUB) << "Allowed Rectangular Region Spanned by:" << endl
        << pad(SUB,w, "  Upper") << upper.x << d << upper.y << d << upper.z << endl
        << pad(SUB,w, "  Lower") << lower.x << d << lower.y << d << lower.z << endl;
      o << indent(SUB) << "Registered Groups:" << endl;
      for (auto g : groups)
        o << indent(SUB) << "  " << g->name << endl;
      return o.str();
    }

    bool RestrictedVolume::outside(const Point &a) {
      if (a.z<lower.z) return true;
      if (a.z>upper.z) return true;
      if (a.x<lower.x) return true;
      if (a.x>upper.x) return true;
      if (a.y<lower.y) return true;
      if (a.y>upper.y) return true;
      return false;
    }

    double RestrictedVolume::g_external(const p_vec &p, Group &g) {
      if (std::find(groups.begin(), groups.end(), &g)!=groups.end())
        for (auto i : g)
          if ( outside( p[i]) )
            return pc::infty;
      return 0;
    }

  }//namespace
}//namespace

