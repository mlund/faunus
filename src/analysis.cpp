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
      stepcnt = 0;
    }

    AnalysisBase::AnalysisBase(Tmjson &j) : w(30), cnt(0), runfraction(1.0) {
      steps = j["steps"] | 1;
      stepcnt = 0;
    }

    AnalysisBase::~AnalysisBase() {}

    /** @todo Remove this function -- let sample() handle it instead */
    bool AnalysisBase::run() {
      if (slump() > runfraction)
        return false;
      cnt++;
      return true;
    }

    void AnalysisBase::_sample() { /* make pure virtual! */ }

    void AnalysisBase::sample() {
      if (steps > stepcnt) {
        timer.start();
        _sample();
        timer.stop();
        stepcnt++;
        cnt++;
      } else stepcnt=0;
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
        o << pad(SUB,w,"Reference:") << cite << "\n";

      o << pad(SUB,w,"Sample interval") << steps << "\n";
      if (steps>1)
        o << pad(SUB,w,"Runfraction") << runfraction*100 << percent << "\n";

      if (cnt>0) {
        o << pad(SUB,w,"Number of sample events") << cnt << "\n";
        double time = timer.result();
        if (time>1e-3)
          o << pad(SUB,w,"Relative time") << time << "\n";
      }
      o << _info();
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
      assert(!g1.name.empty() && "Group 1 (body 1) should have a name.");
      assert(!g2.name.empty() && "Group 2 (body 2) have a name.");
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
      if (f)
        f << p.transpose() << endl;
    }

    string TwobodyForce::_info() {
      using namespace Faunus::textio;
      std::ostringstream o;
      o << "Base class of twobody force." << endl;
      return o.str();
    }

    TwobodyForceDirect::TwobodyForceDirect(InputMap &in, Group &g1, Group &g2, Group &_ions) : TwobodyForce(in, g1, g2, _ions) {
      name="Twobody direct mean force calculation";
      f_pp = f_pi = f_ip = Point(0,0,0);
    }

    Point TwobodyForceDirect::meanforce() {
      Point p = (f_pp+(f_pi-f_ip)*0.5)/(double)cnt;
      return p;
    }

    string TwobodyForceDirect::_info() {
      using namespace Faunus::textio;
      std::ostringstream o;
      o << pad(SUB,w,"Mean direct Force p-p:") << "(" << (f_pp/cnt).transpose() << ") kT/Å" << endl;
      o << pad(SUB,w,"Mean direct Force p-i:") << "(" << (f_pi/cnt).transpose() << ") kT/Å" << endl;
      o << pad(SUB,w,"Mean direct Force i-p:") << "(" << (f_ip/cnt).transpose() << ") kT/Å" << endl;
      o << pad(SUB,w,"Mean direct Force:") << "(" << meanforce().transpose() << ") kT/Å" << endl;
      return o.str();
    }

    TwobodyForceMidp::TwobodyForceMidp(InputMap &in,Group &g1, Group &g2, Group &_ions,
        Analysis::LineDistributionNorm<float,unsigned long int> *_saltdistr) : TwobodyForce(in,g1,g2,_ions) {
      name="Twobody midplane mean force calculation";
      f_pp = f_pi = f_ip = f_ii = Point(0,0,0);
      saltdistr=_saltdistr;
    }

    Point TwobodyForceMidp::meanforce() {
      Point p = f_pp+f_pi-f_ip+f_ii;
      float midp = saltdistr->mid();
      float endp = saltdistr->end();
      Point I = Point(1,1,1);
      return p/(double)cnt+I*(midp-endp);
    }

    string TwobodyForceMidp::_info() {
      using namespace Faunus::textio;
      std::ostringstream o;
      o << pad(SUB,w,"Mean midp Force p-p:") << "(" << (f_pp/cnt).transpose() << ") kT/Å" << endl;
      o << pad(SUB,w,"Mean midp Force p-i:") << "(" << (f_pi/cnt).transpose() << ") kT/Å" << endl;
      o << pad(SUB,w,"Mean midp Force i-p:") << "(" << (f_ip/cnt).transpose() << ") kT/Å" << endl;
      o << pad(SUB,w,"Mean midp Force i-i:") << "(" << (f_ii/cnt).transpose() << ") kT/Å" << endl;
      o << pad(SUB,w,"Midpressure:") << saltdistr->mid() << " kT/Å" << endl;
      o << pad(SUB,w,"Endpressure:") << saltdistr->end() << " kT/Å" << endl;
      o << pad(SUB,w,"Mean midp Force:") << "(" << meanforce().transpose() << ") kT/Å" << endl;
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
      if (!exclusionlist.empty()) {
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
          << setw(k+1) << Z2[m.first].avg()-pow(m.second.avg(),2)
          << setw(k) << mu[m.first].avg()
          << setw(k) << mu2[m.first].avg()-pow(mu[m.first].avg(),2)<< endl;
      return o.str();
    }

    void BilayerStructure::_test(UnitTest &t) {
      t("bilayer_order", S.avg() );
      t("bilayer_area", A.avg() );
    }

    CombinedAnalysis::CombinedAnalysis(AnalysisBase* a, AnalysisBase* b) {
      v.reserve(2);
      v.push_back(a);
      v.push_back(b);
    }

    void CombinedAnalysis::sample() { for (auto i : v) i->sample(); }

    string CombinedAnalysis::info() {
      std::ostringstream o;
      for (auto i : v)
        o << i->info();
      return o.str();
    }

    string CombinedAnalysis::_info() {return string();}

    void CombinedAnalysis::_sample() {}

    CombinedAnalysis& operator+(AnalysisBase &a, AnalysisBase &b) {
      return *(new CombinedAnalysis(&a,&b));
    }
  }//namespace
}//namespace
