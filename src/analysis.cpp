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

    TwobodyForce::TwobodyForce(InputMap &in, Group &g1, Group &g2, Group &_ions) {
      name="Twobody mean force calculation";
      runfraction = in.get<double>("pforce", 1.0);
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

    Point TwobodyForce::meanforce() { return Point(0,0,0); }

    void TwobodyForce::save(string filename) {
      Point p = meanforce();
      std::ofstream f(filename.c_str());
      f.precision(10);
      if (f) {
        f << p.x() << " " << p.y() << " " << p.z() << endl;
      }
    }

    string TwobodyForce::_info() {
      using namespace Faunus::textio;
      std::ostringstream o;
      o << "Base class of twobody force." << endl;
      return o.str();
    }

    TwobodyForceDirect::TwobodyForceDirect(InputMap &in, Group &g1, Group &g2, Group &_ions) : TwobodyForce(in, g1, g2, _ions) {
      name="Twobody direct mean force calculation";
      f_pp = Point(0,0,0);
      f_pi = Point(0,0,0);
      f_ip = Point(0,0,0);
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

    TwobodyForceMidp::TwobodyForceMidp(InputMap &in,Group &g1, Group &g2, Group &_ions,
        Analysis::LineDistributionNorm<float,unsigned long int> *_saltdistr) : TwobodyForce(in,g1,g2,_ions) {
      name="Twobody midplane mean force calculation";
      f_pp = Point(0.0, 0.0, 0.0);
      f_pi = Point(0.0, 0.0, 0.0);
      f_ip = Point(0.0, 0.0, 0.0);
      f_ii = Point(0.0, 0.0, 0.0);
      saltdistr=_saltdistr;
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

    string PolymerShape::_info() {
      char w=10;
      using namespace textio;
      std::ostringstream o;
      o << endl << indent(SUBSUB) << std::left << setw(w) << "Group"
        << setw(w+5) << bracket("Rg"+squared)
        << setw(w+12) << bracket("Rg"+squared)+"-"+bracket("Rg")+squared
        << setw(w+5) << rootof+bracket("Rg"+squared)
        << setw(w+7) << rootof+bracket("Rgx"+squared)
        << setw(w+7) << rootof+bracket("Rgy"+squared)
        << setw(w+7) << rootof+bracket("Rgz"+squared)
        << setw(w+7) << rootof+bracket("Re"+squared)
        << bracket("Re"+squared)+"/"+bracket("Rg"+squared) << endl; 
      for (auto &m : Rg2)
        o << indent(SUBSUB) << std::left << setw(w) << m.first
          << std::setprecision(4)
          << setw(w) << m.second.avg()
          << setw(w+2) << m.second.avg() - pow( Rg[m.first].avg(),2 )
          << setw(w-2) << sqrt( m.second.avg() )
          << setw(w) << sqrt(Rg2x[m.first].avg())
          << setw(w) << sqrt(Rg2y[m.first].avg())
          << setw(w) << sqrt(Rg2z[m.first].avg())
          << setw(w) << sqrt(Re2[m.first].avg())
          << Re2[m.first].avg() / Rg2[m.first].avg() << endl; 
      return o.str();
    }

    void PolymerShape::_test(UnitTest &t) {
      for (auto &m : Rg2)
        t("PolymerShape_Rg"+m.first, Rg[m.first].avg() );
    }

    ChargeMultipole::ChargeMultipole(){
      name="Multipole";
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

    void BilayerStructure::_test(UnitTest &t) {
      t("bilayer_order", S.avg() );
      t("bilayer_area", A.avg() );
    }

  }//namespace
}//namespace
