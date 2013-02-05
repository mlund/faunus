#include <faunus/analysis.h>
#include <faunus/species.h>
#include <faunus/energy.h>
#include <faunus/space.h>
#include <faunus/inputfile.h>
#include <faunus/geometry.h>
#include <faunus/textio.h>

namespace Faunus {

  namespace Analysis {
    AnalysisBase::AnalysisBase() : w(30), cnt(0), runfraction(1.0) {
    }

    AnalysisBase::~AnalysisBase() {}


    bool AnalysisBase::run() {
      if (slp_global() > runfraction)
        return false;
      cnt++;
      return true;
    }

    void AnalysisBase::_test(UnitTest &t) {}

    void AnalysisBase::test(UnitTest &t) {
      _test(t);
    }

    string AnalysisBase::info() {
      assert(!name.empty() && "Please name analysis.");
      using namespace textio;
      std::ostringstream o;
      o << header("Analysis: "+name);
      if (!cite.empty())
        o << pad(SUB,w,"Reference:") << cite << endl;
      o << pad(SUB,w,"Runfraction") << runfraction*100 << percent << endl;
      if (cnt>0)
        o << pad(SUB,w,"Number of sample events") << cnt << endl
          << _info();
      return o.str();
    }
    
    TwobodyForce::TwobodyForce(InputMap &in,Energy::Energybase &e, Space &s, Group &g1, Group &g2, Group &_ions) {
      name="Twobody mean force calculation";
      runfraction = in.get<double>("pforce", 1.0);
      pot=&e;
      spc=&s;
      igroup1=nullptr;
      igroup2=nullptr;
      setTwobodies(g1, g2, _ions);
    }
    
    /*!
     * \brief Set the groups
     * \param g1 Body #1
     * \param g2 Body #2
     * \param _ions Salt, counterions etc that mediates the two bodies
     */
    void TwobodyForce::setTwobodies(Group &g1, Group &g2, Group &_ions) {
      assert(&g1!=nullptr);
      assert(!g1.name.empty() && "Group 1 (body 1) should have a name.");
      assert(&g2!=nullptr);
      assert(!g2.name.empty() && "Group 2 (body 2) have a name.");
      assert(&_ions!=nullptr);
      assert(!_ions.name.empty() && "Group 3 (ions) should have a name.");
      igroup1=&g1;
      igroup2=&g2;
      ions=&_ions;
    }
    
    void TwobodyForce::save(string filename) {
      Point p = meanforce();
      std::ofstream f(filename.c_str());
      f.precision(10);
      if (f) {
        f << p.x() << " " << p.y() << " " << p.z() << endl;
      }
    }
    
    void TwobodyForce::calc() {}
    
    Point TwobodyForce::meanforce() { return Point(0.0, 0.0, 0.0); }
    
    string TwobodyForce::_info() {
      using namespace Faunus::textio;
      std::ostringstream o;
      o << "Base class of twobody force." << endl;
      return o.str();
    }
    
    
    TwobodyForceDirect::TwobodyForceDirect(InputMap &in,Energy::Energybase &e, Space &s, Group &g1, Group &g2, Group &_ions) : TwobodyForce(in, e, s, g1, g2, _ions) {
      name="Twobody direct mean force calculation";
      //f_mean1 = Point(0.0, 0.0, 0.0);
      //f_mean2 = Point(0.0, 0.0, 0.0);
      f_pp = Point(0.0, 0.0, 0.0);
      f_pi = Point(0.0, 0.0, 0.0);
      f_ip = Point(0.0, 0.0, 0.0);
    }
    
    /*!
     * \brief Calculate the direct force between the two bodies
     */
    void TwobodyForceDirect::calc() {
      if (run()) {
        // Force between the two bodies
        for (auto i : *igroup1) {
          for (auto j : *igroup2) {
            Point f = pot->f_p2p(spc->p[i], spc->p[j]);
            f_pp += f;
          }
        }
        //f_pp += 1.0*_f_pp;
        //f_mean1 += 1.0*_f_pp;
        //f_mean2 += -1.0*_f_pp;
        for (auto i : *igroup1) {
          for (auto j : *ions) {
            Point f = pot->f_p2p(spc->p[i], spc->p[j]);
            f_pi += f;
          }
        }
        for (auto i : *igroup2) {
          for (auto j : *ions) {
            Point f = pot->f_p2p(spc->p[i], spc->p[j]);
            f_ip += f;
          }
        }
      }
    }
    
    Point TwobodyForceDirect::meanforce() {
      Point p = (f_pp+(f_pi-f_ip)*0.5)/(double)cnt;
      return p;
    }
    
    string TwobodyForceDirect::_info() {
      using namespace Faunus::textio;
      std::ostringstream o;
      o << pad(SUB,w,"Mean direct Force p-p:") << "(" << f_pp.x()/(double)cnt << ", " << f_pp.y()/(double)cnt << ", " << f_pp.z()/(double)cnt << ") kT/Å" << endl;
      o << pad(SUB,w,"Mean direct Force p-i:") << "(" << f_pi.x()/(double)cnt << ", " << f_pi.y()/(double)cnt << ", " << f_pi.z()/(double)cnt << ") kT/Å" << endl;
      o << pad(SUB,w,"Mean direct Force i-p:") << "(" << f_ip.x()/(double)cnt << ", " << f_ip.y()/(double)cnt << ", " << f_ip.z()/(double)cnt << ") kT/Å" << endl;
      o << pad(SUB,w,"Mean direct Force:") << "(" << meanforce().x() << ", " << meanforce().y() << ", " << meanforce().z() << ") kT/Å" << endl;
      return o.str();
    }
    
    
    TwobodyForceMidp::TwobodyForceMidp(InputMap &in,Energy::Energybase &e, Space &s, Group &g1, Group &g2, Group &_ions, Analysis::LineDistributionNorm<float,unsigned long int> *_saltdistr) : TwobodyForce(in, e, s, g1, g2, _ions) {
      name="Twobody midplane mean force calculation";
      f_pp = Point(0.0, 0.0, 0.0);
      f_pi = Point(0.0, 0.0, 0.0);
      f_ip = Point(0.0, 0.0, 0.0);
      f_ii = Point(0.0, 0.0, 0.0);
      saltdistr=_saltdistr;
    }
    
    /*!
     * \brief Calculate the direct force between the two bodies
     */
    void TwobodyForceMidp::calc() {
      if (run()) {
        // Force between the two bodies
        for (auto i : *igroup1) {
          for (auto j : *igroup2) {
            Point f = pot->f_p2p(spc->p[i], spc->p[j]);
            f_pp += f;
          }
        }
        
        for (auto i : *igroup1) {
          for (auto j : *ions) {
            if (spc->p[j].z() < 0.0) {
              Point f = pot->f_p2p(spc->p[i], spc->p[j]);
              f_pi += f;
            }
          }
        }
        
        for (auto i : *igroup2) {
          for (auto j : *ions) {
            if (spc->p[j].z() >= 0.0) {
              Point f = pot->f_p2p(spc->p[i], spc->p[j]);
              f_ip += f;
            }
          }
        }
        
        for (auto i : *ions) {
          if (spc->p[i].z() >= 0.0) {
            for (auto j : *ions) {
              if (spc->p[j].z() < 0.0) {
                Point f = pot->f_p2p(spc->p[i], spc->p[j]);
                f_ii += f;
              }
            }
          }
        }
      }
    }
    
    Point TwobodyForceMidp::meanforce() {
      Point p = f_pp+f_pi-f_ip+f_ii;
      float midp = saltdistr->mid();
      float endp = saltdistr->end();
      Point I = Point(1.0, 1.0, 1.0);
      return p/(double)cnt+I*(midp-endp);
    }
    
    string TwobodyForceMidp::_info() {
      using namespace Faunus::textio;
      std::ostringstream o;
      o << pad(SUB,w,"Mean midp Force p-p:") << "(" << f_pp.x()/(double)cnt << ", " << f_pp.y()/(double)cnt << ", " << f_pp.z()/(double)cnt << ") kT/Å" << endl;
      o << pad(SUB,w,"Mean midp Force p-i:") << "(" << f_pi.x()/(double)cnt << ", " << f_pi.y()/(double)cnt << ", " << f_pi.z()/(double)cnt << ") kT/Å" << endl;
      o << pad(SUB,w,"Mean midp Force i-p:") << "(" << f_ip.x()/(double)cnt << ", " << f_ip.y()/(double)cnt << ", " << f_ip.z()/(double)cnt << ") kT/Å" << endl;
      o << pad(SUB,w,"Mean midp Force i-i:") << "(" << f_ii.x()/(double)cnt << ", " << f_ii.y()/(double)cnt << ", " << f_ii.z()/(double)cnt << ") kT/Å" << endl;
      o << pad(SUB,w,"Midpressure:") << saltdistr->mid() << " kT/Å" << endl;
      o << pad(SUB,w,"Endpressure:") << saltdistr->end() << " kT/Å" << endl;
      o << pad(SUB,w,"Mean midp Force:") << "(" << meanforce().x() << ", " << meanforce().y() << ", " << meanforce().z() << ") kT/Å" << endl;
      return o.str();
    }
    
    
    PolymerShape::PolymerShape() {
      name="Polymer Shape";
    }

    Point PolymerShape::vectorgyrationRadiusSquared(const Group &pol, const Space &spc) {
      assert( spc.geo->dist(pol.cm, pol.massCenter(spc))<1e-9 && "Mass center must be in sync.");
      double sum=0;
      Point t, r2(0,0,0);
      for (auto i : pol) {
        t = spc.p[i]-pol.cm;                // vector to center of mass
        spc.geo->boundary(t);               // periodic boundary (if any)
        r2.x() += spc.p[i].mw * t.x() * t.x();
        r2.y() += spc.p[i].mw * t.y() * t.y();
        r2.z() += spc.p[i].mw * t.z() * t.z();
        sum += spc.p[i].mw;                 // total mass
      }
      assert(sum>0 && "Zero molecular weight not allowed.");
      return r2*(1./sum);
    }

    double PolymerShape::gyrationRadiusSquared(const Group &pol, const Space &spc) {
      assert( spc.geo->dist(pol.cm, pol.massCenter(spc))<1e-9 && "Mass center must be in sync.");
      Point rg2=vectorgyrationRadiusSquared(pol,spc);
      return rg2.x()+rg2.y()+rg2.z();
    }

    Point PolymerShape::vectorEnd2end(const Group &pol, const Space &spc) {
      return spc.geo->vdist( spc.p[pol.front()], spc.p[pol.back()] );
    }

    void PolymerShape::sample(const Group &pol, const Space &spc) {
      if (!run() || pol.front()==pol.back())
        return;
      Point r2 = vectorgyrationRadiusSquared(pol,spc);
      double rg2 = r2.x()+r2.y()+r2.z(); 
      double re2 = spc.geo->sqdist( spc.p[pol.front()], spc.p[pol.back()] );
      Rg2[pol.name]  += rg2;
      Rg2x[pol.name] += r2.x();
      Rg2y[pol.name] += r2.y();
      Rg2z[pol.name] += r2.z();
      Rg[pol.name]   += sqrt(r2.x()+r2.y()+r2.z());
      Re2[pol.name]  += re2; //end-2-end squared
      double rs = Re2[pol.name].avg()/Rg2[pol.name].avg(); // fluctuations in shape factor
      Rs[pol.name]   += rs;
      Rs2[pol.name]  += rs*rs;
      
      //Point re = vectorEnd2end(pol,spc);
      //Re2[pol.name] += pow(re.len(), 2);
    }

    string PolymerShape::_info() {
      char w=13;
      using namespace textio;
      std::ostringstream o;
      o << endl << indent(SUBSUB) << std::left << setw(w) << "Polymer"
        << setw(w+5) << bracket("Rg"+squared)
        << setw(w+12) << bracket("Rg"+squared)+"-"+bracket("Rg")+squared
        << setw(w+7) << rootof+bracket("Rg"+squared)
        << setw(w+7) << rootof+bracket("Rgx"+squared)
        << setw(w+6) << rootof+bracket("Rgy"+squared)
        << setw(w+6) << rootof+bracket("Rgz"+squared)
        << setw(w+7) << rootof+bracket("Re"+squared)
        << setw(w+7) << bracket("Re"+squared)+"/"+bracket("Rg"+squared) << endl; 
      for (auto &m : Rg2)
        o << indent(SUBSUB) << std::left << setw(w) << m.first << setw(w) << m.second.avg()
          << setw(w) << m.second.avg() - pow( Rg[m.first].avg(),2 )
          << setw(w) << sqrt( m.second.avg() )
          << setw(w) << sqrt(Rg2x[m.first].avg())
          << setw(w) << sqrt(Rg2y[m.first].avg())
          << setw(w) << sqrt(Rg2z[m.first].avg())
          << setw(w) << sqrt(Re2[m.first].avg())
          << setw(w) << Re2[m.first].avg() / Rg2[m.first].avg() << endl; 
      return o.str();
    }

    void PolymerShape::_test(UnitTest &t) {
      for (auto &m : Rg2)
        t("PolymerShape_Rg"+m.first, Rg[m.first].avg() );
    }

    ChargeMultipole::ChargeMultipole(){
      name="Multipole";
    }

    /*!
     * \param g Group to calculate charge for
       \param spc Space
     */
    double ChargeMultipole::charge(const Group &g,const Space &spc) {
      double x=0;
      for (auto i : g){
        if (exclude(spc.p[i])==false)
          x+=spc.p[i].charge;
      }
      return x;
    }

    double ChargeMultipole::dipole(const Group &g, const Space &spc){
      assert( spc.geo->dist(g.cm, g.massCenter(spc))<1e-9 && "Mass center must be in sync.");
      Point t, mu(0,0,0);
      for (auto i : g) {
        if (exclude(spc.p[i])==false){
          t = spc.p[i]-g.cm;                // vector to center of mass
          spc.geo->boundary(t);               // periodic boundary (if any)
          mu.x()+=spc.p[i].charge * t.x();
          mu.y()+=spc.p[i].charge * t.y();
          mu.z()+=spc.p[i].charge * t.z();
        }
      }
      return mu.len();
    }

    bool ChargeMultipole::exclude(const particle &p){
      if (exclusionlist.find(atom[p.id].name)==exclusionlist.end())
        return false;
      return true;
    }

    void ChargeMultipole::sample(const Group &g, const Space &spc) {
      assert(!g.name.empty() && "All Groups should have a name!");
      if (!run())
        return;
      double z=charge(g, spc);
      Z[g.name]+=z;
      Z2[g.name]+=pow(z,2);
      double dip=dipole(g,spc);
      mu[g.name]+=dip;
      mu2[g.name]+=pow(dip,2);
    }

    string ChargeMultipole::_info(){
      using namespace textio;
      char k=13;
      std::ostringstream o;
      if (~exclusionlist.empty()) {
        o << pad(SUB,w, "Exclusion list");
        for (auto i : exclusionlist)
          o << i << " ";
      }
      o << endl << indent(SUB) << std::left << setw(w) << "Macromolecule  "
        << setw(k+4) << bracket("Z")
        << setw(k+11) << bracket("Z"+squared)+"-"+bracket("Z")+squared
        << setw(k+5) << bracket(textio::mu)
        << setw(k+5) << bracket(textio::mu+squared)+"-"+bracket(textio::mu)+squared << endl;
      for (auto &m : Z)
        o << indent(SUB) << std::left << setw(w) << m.first << setw(k) << m.second.avg()
          << setw(k) << Z2[m.first].avg()-pow(m.second.avg(),2)
          << setw(k) << mu[m.first].avg()
          << setw(k) << mu2[m.first].avg()-pow(mu[m.first].avg(),2)<< endl;
      return o.str();
    }

    Widom::Widom(Space &spc, Energy::Energybase &pot) {
      name="Multi Particle Widom Analysis";
      cite="doi:10/dkv4s6";
      spcPtr=&spc;
      potPtr=&pot;
    }

    void Widom::sample(int ghostin) {
      if (!run())
        return;
      assert(spcPtr->geo!=NULL);
      int n=g.size();
      for (int k=0; k<ghostin; k++) {     // insert ghostin times
        double du=0;
        for (int i=0; i<n; i++)
          spcPtr->geo->randompos( g[i] ); // random ghost positions
        for (int i=0; i<n; i++)
          du+=potPtr->all2p( spcPtr->p, g[i] );    // energy with all particles in space
        for (int i=0; i<n-1; i++)
          for (int j=i+1; j<n; j++)
            du+=potPtr->p2p( g[i], g[j] );   // energy between ghost particles
        expsum += exp(-du);
      }
    }

    void Widom::addGhost(particle p) {
      g.push_back(p);
    }

    void Widom::addGhost(Space &c) {
      std::map<short,bool> map;
      for (auto p : c.p)
        map[ p.id ] = true;
      for (auto &m : map) {
        particle a;
        a=atom[m.first];
        addGhost(a);
      }
    }

    void Widom::check(UnitTest &test) {
      test("widom_muex", muex() );
    }

    string Widom::_info() {
      using namespace Faunus::textio;
      std::ostringstream o;
      o << pad(SUB,w, "Number of insertions") << expsum.cnt << endl
        << pad(SUB,w, "Excess chemical pot.") << muex() << kT << endl
        << pad(SUB,w, "Mean activity coefficient") << gamma() << endl
        << pad(SUB,w, "Ghost particles");
      for (auto p : g)
        o << atom[p.id].name << " ";
      o << endl;
      return o.str();
    }

    double Widom::gamma() {
      return exp(muex());
    }

    double Widom::muex() {
      return -log(expsum.avg())/g.size();
    }

    //--------------------------------------------------------------------------------

    /**
     * @param bjerrumLength Bjerrum length [angstrom]
     * @param insertions Number of insertions per call to `insert()`
     */
    WidomScaled::WidomScaled(double bjerrumLength, int insertions) {
      assert(insertions>=0);
      assert(bjerrumLength>0);
      name="Single particle Widom insertion w. charge scaling"; 
      cite="doi:10/ft9bv9 + doi:10/dkv4s6"; 
      lB=bjerrumLength;
      ghostin=insertions;
    }

    /**
     * This will add particle `p` to the list of ghost particles
     * to insert.
     */
    void WidomScaled::add(const particle &p) {
      g.push_back(p);
      init();
    }

    /**
     * This will scan the particle vector for particles and each unique type
     * will be added to the list a ghost particles to insert.
     */
    void WidomScaled::add(const p_vec &p) {
      std::set<particle::Tid> ids;
      for (auto &i : p)
        ids.insert(i.id);
      for (auto i : ids) {
        particle a;
        a=atom[i];
        add(a);
      }
    }

    bool WidomScaled::overlap(const particle &a, const particle &b, const Geometry::Geometrybase &geo)
    {
      double s=a.radius+b.radius;
      return (geo.sqdist(a,b)<s*s) ? true : false;
    }

    void WidomScaled::init() {
      int gspec=g.size();
      chel.resize(gspec);
      chhc.resize(gspec);
      chex.resize(gspec);
      chtot.resize(gspec);
      ewden.resize(gspec);
      ewnom.resize(gspec);
      chint.resize(gspec);
      expuw.resize(gspec);
      chexw.resize(gspec);
      ihc.resize(gspec);
      irej.resize(gspec);

      for (int i=0; i<gspec; i++){
        chel[i]=0;
        chhc[i]=0;
        chex[i]=0;
        chtot[i]=0;
        ihc[i]=0;
        ewden[i].resize(11);
        ewnom[i].resize(11);
        chint[i].resize(11);
        for(int j=0; j<11; j++)
          ewden[i][j] = ewnom[i][j] = chint[i][j] = 0;
      }
    }

    /**
     * @param p List of particles to insert into. This will typically be the main
     *          particle vector, i.e. `Space::p`.
     * @param geo Geometry to use for distance calculations and random position generation
     */
    void WidomScaled::insert(const p_vec &p, Geometry::Geometrybase &geo) {
      assert(lB>0);
      assert(&geo!=nullptr);
      if (!run() || g.empty() || p.empty())
        return;
      particle ghost;
      double u,cu;
      for (int i=0; i<ghostin; i++) {
        geo.randompos(ghost);
        int goverlap=0;
        for (size_t k=0; k<g.size(); k++) {
          ghost.radius = g[k].radius;
          irej[k]=0;
          int j=0;
          while (!overlap(ghost,p[j],geo) && j<(int)p.size())
            j++;
          if (j!=(int)p.size()) {
            ihc[k]++;
            irej[k]=1;
            goverlap++;
          }
        }

        if ( goverlap != (int)g.size() ) {
          cu=0;
          u=0;  //elelectric potential (Coulomb only!)
          for (auto &i : p) {
            double invdi=1/geo.dist(ghost,i);
            cu+=invdi;
            u+=invdi*i.charge;
          } 
          cu=cu*lB;
          u=u*lB;
          double ew,ewla,ewd;
          for (size_t k=0; k < g.size(); k++) {
            if (irej[k]==0) {
              expuw[k]+=exp(-u*g[k].charge);
              for (int cint=0; cint<11; cint++) {
                ew=g[k].charge*(u-double(cint)*0.1*g[k].charge*cu/double(p.size()));
                ewla = ew*double(cint)*0.1;
                ewd=exp(-ewla);
                ewden[k][cint]+=ewd;
                ewnom[k][cint]+=ew*ewd;
              }
            }
          }
        }
      }
    }

    string WidomScaled::_info() {
      using namespace textio;
      std::ostringstream o;
      double aint4,aint2,aint1;
      for(size_t i=0; i<g.size(); i++) {
        for(int cint=0; cint<11; cint++) {
          if(ewden[i][cint]==0)
            std::cout << "# WARNING: Widom denominator equals zero" << endl;
          else
            chint[i][cint]=ewnom[i][cint]/ewden[i][cint];
        }
        aint4=chint[i][1]+chint[i][3]+chint[i][5]+chint[i][7]+chint[i][9];
        aint2=chint[i][2]+chint[i][4]+chint[i][6]+chint[i][8];
        aint1=chint[i][0]+chint[i][10];  
        chel[i]=1./30.*(aint1+2*aint2+4*aint4);
      }

      double cnttot;
      cnttot=cnt*ghostin;
      o << pad(SUB,w,"Number of Insertions") << cnttot << endl
        << pad(SUB,w,"Excess chemical potentials (kT)") << endl
        << "            total   elec.  hs            z       r"<< endl;
      for (size_t i=0; i < g.size(); i++) {
        chhc[i]=-log(double(cnttot-ihc[i])/cnttot);
        chexw[i]=-log(expuw[i]);
        chex[i]=chhc[i]+chel[i];
        o.unsetf( std::ios_base::floatfield );
        o << "    [" << i << "] "
          << std::setprecision(4)
          << std::setw(9) << chex[i]
          << std::setw(9) << chel[i]
          << std::setw(9) << chhc[i]
          << std::setprecision(2) << std::fixed
          << std::setw(9) << g[i].charge
          << std::setw(9) << g[i].radius << endl;
      }
      return o.str();
    }
  }//namespace
}//namespace
