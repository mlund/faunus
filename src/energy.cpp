#include <faunus/common.h>
#include <faunus/inputfile.h>
#include <faunus/point.h>
#include <faunus/energy.h>
#include <faunus/textio.h>
#include <faunus/space.h>
#include <faunus/species.h>

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

    /*!
     * This function sets the volume of the Geometry used for the energy calculations.
     * Currently this is accessed through the geo pointer, but in the future each Energybase
     * derivative who needs a Geometry should have it's own instance of a Geometry.
     */
    void Energybase::setVolume(double vol) {
      if (geo!=nullptr)
        geo->setVolume(vol);
    }

    void Energybase::setTemperature(double) {
      assert(!"Not yet implemented!");
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

    ExternalPressure::ExternalPressure(Geometry::Geometrybase &e, double pressure) {
      assert(&e!=nullptr && "Geometry must be defined for this energy term");
      assert(pressure>=0 && "Specify non-negative pressure");
      name="External Pressure";
      setGeometry(e);
      P=pressure;
    }

    double ExternalPressure::external() {
      double V=geo->getVolume();
      assert(V>1e-9 && "Volume must be non-zero!");
      return P*V - log(V); 
    }

    double ExternalPressure::g_external(const p_vec &p, Group &g) {
      double N=1,
             V=geo->getVolume();
      if (g.isAtomic())
        N=g.size();
      return -N*log(V);
    }

    string ExternalPressure::_info() {
      char w=15;
      using namespace textio;
      std::ostringstream o;
      o << pad(SUB,w,"Pressure") << P*1e30/pc::Nav << " mM = " << P*pc::kB*pc::T()*1e30 << " Pa" << endl;
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
          assert(!"error");
          std::cerr << s;
          return;
        }
        e.setGeometry(*geo);
      }
      else if (geo==nullptr)
        geo = &e.getGeometry();
      baselist.push_back( &e );
    }

    void Hamiltonian::setVolume(double vol) {
      for (auto e : baselist )
        e->setVolume(vol);
    }

    void Hamiltonian::setTemperature(double T) {
      for (auto e : baselist )
        e->setTemperature(T);
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

    /*!
     * For early rejection we do a reverse loop over energy
     * classes to test if the energy is infinity. Reverse
     * because constraining energy classes are often added
     * last.
     */
    double Hamiltonian::g_external(const p_vec &p, Group &g) {
      double u=0;
      for (auto b=baselist.rbegin(); b!=baselist.rend(); ++b) {
        u += (*b)->g_external(p,g);
        if (u>=pc::infty)
          break;
      }
      return u;
      //for (auto b : baselist)
      //  u += b->g_external(p,g);
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

    Bonded::Bonded() {
      name="Bonded particles";
      geo=nullptr;
    }

    Bonded::Bonded(Geometry::Geometrybase &g) {
      name="Bonded particles";
      geo=&g;
    }

    string Bonded::_info() {
      using namespace Faunus::textio;
      std::ostringstream o;
      o << indent(SUBSUB) << std::left
        << setw(7) << "i" << setw(7) << "j" << endl;
      for (auto &m : list)
        o << indent(SUBSUB) << std::left << setw(7) << m.first.first
          << setw(7) << m.first.second << m.second->brief() << endl;
      return o.str();
    }

    /*!
     * \todo Optimize by using iterator directly as returned by find.
     */
    double Bonded::i2i(const p_vec &p, int i, int j) {
      assert(i!=j);
      auto f=list.find( pair_permutable<int>(i,j) );
      if (f!=list.end())
        return f->second->tokT() * f->second->operator()( p[i], p[j], geo->sqdist( p[i], p[j] ) );
      return 0;
    }

    double Bonded::i2all(const p_vec &p, int i) {
      assert(geo!=nullptr);  //debug
      assert( i>=0 && i<(int)p.size() ); //debug

      double u=0;
      for (auto &m : list) {
        int j=m.first.first;
        int k=m.first.second;
        assert(j!=k && "Pairs between identical atom index not allowed.");
        if (i==j || i==k)
          u+=m.second->tokT() * m.second->operator()( p[j], p[k], geo->sqdist( p[j], p[k] ) );
      }
      return u;
    }

    double Bonded::total(const p_vec &p) {
      assert(&geo!=nullptr);  //debug
      double u=0;
      for (auto &m : list) {
        int i=m.first.first;
        int j=m.first.second;
        assert(i!=j && "Pairs between identical atom index not allowed.");
        assert(i>=0 && i<(int)p.size() && j>=0 && j<(int)p.size()); //debug
        u += m.second->tokT() * m.second->operator()( p[i], p[j], geo->sqdist( p[i], p[j] ) );
      }
      return u;
    }

    double Bonded::g2g(const p_vec &p, Group &g1, Group &g2) {
      double u=0;
      for (auto &m : list) {
        int i=m.first.first;
        int j=m.first.second;
        assert(i!=j && "Pairs between identical atom index not allowed.");
        if (g1.find(i))
          if (g2.find(j))
            u+=m.second->tokT() * m.second->operator()( p[i], p[j], geo->sqdist( p[i], p[j] ) );
        if (g1.find(j))
          if (g2.find(i))
            u+=m.second->tokT() * m.second->operator()( p[i], p[j], geo->sqdist( p[i], p[j] ) );
      }
      return u;
    }

    /*!
     * Accounts for bonds between particles within a group. Bonds with particles
     * outside the group are skipped and should be accounted for by the g2g() energy function.
     */
    double Bonded::g_internal(const p_vec &p, Group &g) {
      assert(geo!=nullptr);  //debug
      double u=0;
      for (auto &m : list) {
        int i=m.first.first;
        int j=m.first.second;
        assert(i!=j && "Pairs between identical atom index not allowed.");
        assert(i>=0 && i<(int)p.size() && j>=0 && j<(int)p.size()); //debug
        if (g.find(i))
          if (g.find(j))
            u += m.second->tokT() * m.second->operator()( p[i],p[j],geo->sqdist( p[i],p[j] ) );
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
    RestrictedVolume::RestrictedVolume(InputMap &in, string prefix) {
      name="Restricted Volume";
      lower.x = in.get<double>(prefix+"_lower.x", -pc::infty);
      lower.y = in.get<double>(prefix+"_lower.y", -pc::infty);
      lower.z = in.get<double>(prefix+"_lower.z", -pc::infty);
      upper.x = in.get<double>(prefix+"_upper.x", pc::infty);
      upper.y = in.get<double>(prefix+"_upper.y", pc::infty);
      upper.z = in.get<double>(prefix+"_upper.z", pc::infty);
      assert(upper.x>lower.x && "Upper bound must be bigger than lower bound!");
      assert(upper.y>lower.y && "Upper bound must be bigger than lower bound!");
      assert(upper.z>lower.z && "Upper bound must be bigger than lower bound!");
    }

    string RestrictedVolume::_info() {
      using namespace textio;
      char w=20;
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
      if (a.x<lower.x) return true;
      if (a.y<lower.y) return true;
      if (a.z<lower.z) return true;
      if (a.x>upper.x) return true;
      if (a.y>upper.y) return true;
      if (a.z>upper.z) return true;
      return false;
    }

    double RestrictedVolume::g_external(const p_vec &p, Group &g) {
      if (std::find(groups.begin(), groups.end(), &g)!=groups.end())
        for (auto i : g)
          if ( outside( p[i]) )
            return pc::infty;
      return 0;
    }

    RestrictedVolumeCM::RestrictedVolumeCM(InputMap &in, string prefix) : RestrictedVolume(in,prefix) {
      name+=" (mass center)";
    }

    /*!
     * This will calculate the mass center of the groups based on the
     * given particle vector. If Outside the allowed region, infinity
     * is returned - zero otherwise.
     * \note The Group object is not changed - i.e. cm/cm_trial remain
     * intact.
     */
    double RestrictedVolumeCM::g_external(const p_vec &p, Group &g) {
      if (std::find(groups.begin(), groups.end(), &g)!=groups.end())
        if ( outside( Geometry::massCenter(*geo, p, g) ) )
          return pc::infty;
      return 0;
    }

    /*!
     * If a single particle is moved, this function will check if
     * that particle is within a constrained groups. If so, the
     * above g_external() is called.
     */
    double RestrictedVolumeCM::i_external(const p_vec &p, int i) {
      double u=0;
      for (auto g : groups)     // loop over all groups in space
        if (g->find(i))         // does i belong to one of them?
          u+=g_external(p, *g); // if so, check mass center constraint
      return u;
    }

    void MassCenterConstrain::addPair(Group &a, Group &b, double mindist, double maxdist) {
      data d = {mindist, maxdist};
      pair_permutable<Group*> p(&a, &b);
      gmap[p] = d;
    }

    MassCenterConstrain::MassCenterConstrain(Geometry::Geometrybase &geo) {
      name="Group Mass Center Distance Constraints";
      setGeometry(geo);
    }

    double MassCenterConstrain::g_external(const p_vec &p, Group &g1) {
      for (auto m : gmap)              // scan through pair map
        if (m.first.find(&g1)) {       // and look for group g1
          Group *g2ptr;                // pointer to constrained partner
          if (&g1 == m.first.first)
            g2ptr = m.first.second;
          else
            g2ptr = m.first.first;
          Point cma = Geometry::massCenter(*geo, p, g1);
          Point cmb = Geometry::massCenter(*geo, p, *g2ptr);
          double r2 = geo->sqdist(cma,cmb);
          double min = m.second.mindist;
          double max = m.second.maxdist;
          if (r2<min*min || r2>max*max) 
            return pc::infty;
        }
      return 0;
    }

    string MassCenterConstrain::_info() {
      using namespace Faunus::textio;
      std::ostringstream o;
      o << indent(SUB) << "The following groups have mass center constraints:\n";
      for (auto m : gmap)
        o << indent(SUBSUB) << m.first.first->name << " " << m.first.second->name
          << " " << m.second.mindist << "-" << m.second.maxdist << _angstrom << endl;
      return o.str();
    }

    /*!
     * The InputMap is searched for the following keywords:
     * \li \c dh_ionicstrength - via Potential::DebyeHuckel
     * \li \c gouychapman_phi0 - surface potential [V]
     * \li \c gouychapman_qarea - surface charge density (if phi0 not defined)
     * \li \c gouychapman_rho - surface charge density [1/A^2] (if qarea not defined)
     * \li \c gouychapman_linearize - set to yes for linearized PB (default: no)
     *
     * Equations:
     * \f[ \rho = \sqrt\frac{2 c_0}{\pi l_B}  \sinh( \beta \phi_0 e / 2 ) \f]
     * \f[ \beta \phi_0 e = 2\mbox{~asinh} \left ( \rho \sqrt\frac{\pi l_B} {2 c_0} \right ) \f]
     * \f[ \Gamma_0=\tanh{ \beta \phi_0 z e / 4 } \f]
     * where \f$ lB \f$ is the Bjerrum length, \f$\kappa\f$ the inverse Debye length, and \f$c_0\f$ the
     * bulk salt concentration.
     */
    GouyChapman::GouyChapman(InputMap &in) : dh(in) {
      name = "Gouy-Chapman External Potential";
      string prefix="gouychapman_";
      zposPtr=nullptr;
      c0=dh.ionicStrength() * pc::Nav / 1e27;     // assuming 1:1 salt, so c0=I
      lB=dh.bjerrumLength();
      kappa=1/dh.debyeLength();

      linearize = in.get<bool>(prefix+"linearize",false);

      phi0=in.get<double>(prefix+"phi0",0);     // Surface potential [unitless] = beta * e * phi0
      if ( std::abs(phi0)>1e-6 ) {
        rho=sqrt(2*c0/(pc::pi*lB))*sinh(.5*phi0);   // [Evans & Wennerström, 1999, Colloidal Domain p 138-140]
      }
      else {
        rho=1/in.get<double>(prefix+"qarea",0);
        if (rho>1e20)
          rho=in.get<double>(prefix+"rho",0);
        phi0=2.*asinh(rho * sqrt(.5*lB*pc::pi/c0 ));//  [Evans..]
      }

      gamma0=tanh(phi0/4);                          // assuming z=1  [Evans..]
    }

    string GouyChapman::_info() {
      using namespace textio;
      char w=30;
      std::ostringstream o;
      o << pad(SUB,w,"Surface z-position") << *zposPtr << _angstrom << endl
        << pad(SUB,w,"Bjerrum length") << lB << _angstrom << endl
        << pad(SUB,w,"Debye length") << 1./kappa << _angstrom  << endl
        << pad(SUB,w,"Ionic strenght") << dh.ionicStrength()*1e3<< " mM" << endl
        << pad(SUB,w,"Bulk 1:1 salt concentration") << c0 << _angstrom+cubed << endl
        << pad(SUB,w,"Surface potential") << phi0*pc::kB*pc::T()/pc::e << " J/C=volts" << endl
        << pad(SUB,w,"Unitless surface potential") << phi0 << endl
        << pad(SUB,w,"Area per surface charge") << 1/rho << _angstrom+squared << endl
        << pad(SUB,w,"Surface charge density") << rho*pc::e*1e20  << " C/m" + squared << endl
        << pad(SUB,w,"GC-coefficient Gamma_0") << gamma0  << "  " << endl
        << pad(SUB,w,"Linearized PB") << ((linearize) ? "yes" : "no") << endl;
      return o.str();
    }

    /*!
     * Here, give a reference to the z-position of the Gouy-Chapman surface. A
     * typical usage would be to use len_half.z from Geometry::Cuboidslit which will
     * place the surface at the edge of the cuboid container. By using a reference
     * (pointer) we ensure that the position will follow possible volume fluctuations
     * when in the NPT ensemble (although Poisson-Boltzmann electrostatics may require
     * further considerations when using constant pressure).
     */
    void GouyChapman::setPosition(double &zposition) {
      zposPtr=&zposition;
    }

    double GouyChapman::i_external(const p_vec &p, int i) {
      return p_external(p[i]);
    }

    double GouyChapman::g_external(const p_vec &p, Group &g) {
      double u=0;
      for (auto i : g)
        u+=p_external(p[i]);
      return u;
    }

    MeanFieldCorrection::MeanFieldCorrection(InputMap& in)
      : bin(in.get<double>("mfc_binsize", 2)), 
        dh(in),
        qdensity(bin,Ttable::XYDATA) {
      name = "Mean Field Correction";  
      threshold=in.get<double>("cylinder_radius",pc::infty);
      loadfromdisk=in.get<bool>("mfc_load", false);
      filename=textio::prefix+"mfc_qdensity";
      prefactor=std::exp(-threshold/dh.debyeLength())*dh.bjerrumLength()*pc::pi*2*bin*dh.debyeLength();
      if (loadfromdisk)
        qdensity.load(filename);
    }

    double MeanFieldCorrection::i_external(const p_vec& p, int i) {
      return p[i].charge*qdensity(p[i].z)*prefactor;
    }

    double MeanFieldCorrection::g_external(const p_vec& p, Group& g) {
      double u=0;
      for (auto i : g)  //i will be the index of the group
        u+=i_external(p, i);
      return u;
    }

    /*!
     * Sampling is done only when \c loadfromdisk is set to false. After each sampling event, the charge
     * density table is saved to disk.
     */
    void MeanFieldCorrection::sample(const p_vec& p, double zmin, double zmax){
      if (!loadfromdisk) {
        typedef Analysis::Table2D<double,double> T;
        T qsum(bin,T::XYDATA);  // summed charge in each bin
        double dV=pc::pi*pow(threshold, 2)*bin; // volume of each bin
        for (auto &i : p)
          qsum(i.z)+=i.charge;
        for (double z=zmin; z<=zmax; z+=bin)
          qdensity(z)+=qsum(z)/dV;
        qdensity.save(filename);
      }
    }

    string MeanFieldCorrection::_info(){
      using namespace textio;
      char w=30;
      std::ostringstream o;
      o << pad(SUB,w,"Mean Field hole radius") << threshold << _angstrom << endl
        << pad(SUB,w,"Mean Field bin width") << bin << _angstrom << endl
        << pad(SUB,w,"Prefactor") << prefactor << _angstrom+cubed << endl;
      return o.str();
    }
    
    /*
       pair_permutable<particle::Thydrophobic> createPairHydrophobic(const particle &a, const particle &b) {
       return pair_permutable<particle::Thydrophobic>(a.hydrophobic, b.hydrophobic);
       }
       */

    PairListID::PairListID() : GeneralPairList<particle::Tid>(makepair) {
      name+=" (particle id's)";
    }

    pair_permutable<particle::Tid> PairListID::makepair(const particle &a, const particle &b) {
      return pair_permutable<particle::Tid>(a.id, b.id);
    }

    string PairListID::_info() {
      using namespace Faunus::textio;
      std::ostringstream o;
      o << indent(SUBSUB) << std::left
        << setw(7) << "i" << setw(7) << "j" << endl;
      for (auto &m : list)
        o << indent(SUBSUB) << std::left << setw(7) << atom[m.first.first].name
          << setw(7) << atom[m.first.second].name << m.second->brief() << endl;
      return o.str();
    }

    PairListHydrophobic::PairListHydrophobic() : GeneralPairList<particle::Thydrophobic>(makepair) {
      name+=" (hydrophobic particles)";
    }

    pair_permutable<particle::Thydrophobic> PairListHydrophobic::makepair(const particle &a, const particle &b) {
      return pair_permutable<particle::Thydrophobic>(a.hydrophobic, b.hydrophobic);
    }

    string PairListHydrophobic::_info() {
      using namespace Faunus::textio;
      std::ostringstream o;
      o << indent(SUBSUB) << std::left
        << setw(7) << "i" << setw(7) << "j" << endl;
      for (auto &m : list)
        o << indent(SUBSUB) << std::left << setw(7) << ((m.first.first==true) ? "true":"false")
          << setw(7) << ((m.first.second==true) ? "true":"false") << m.second->brief() << endl;
      return o.str();
    }

    DebyeHuckelActivity::DebyeHuckelActivity(InputMap &in) : dh(in) {
      name="Debye-Huckel Activity Coefficients";
    }

    string DebyeHuckelActivity::_info() { return dh.brief(); }

    double DebyeHuckelActivity::p_external(const particle &p) {
      return dh.excessChemPot(p.charge, 2*p.radius);
    }

    double DebyeHuckelActivity::i_external(const p_vec &p, int i) {
      return p_external(p[i]);
    }

    double DebyeHuckelActivity::g_external(const p_vec &p, Group &g) {
      double u=0;
      for (auto i : g)
        u+=p_external(p[i]);
      return u;
    }

    double systemEnergy(Space &spc, Energy::Energybase &pot, const p_vec &p) {
      double u = pot.external();
      for (auto g : spc.groupList())
        u += pot.g_external(p, *g) + pot.g_internal(p, *g);
      for (size_t i=0; i<spc.groupList().size()-1; i++)
        for (size_t j=i+1; j<spc.groupList().size(); j++)
          u += pot.g2g(p, *spc.groupList()[i], *spc.groupList()[j]);
      return u;
    }

  }//namespace
}//namespace

