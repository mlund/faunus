#include <faunus/analysis.h>
#include <faunus/species.h>
#include <faunus/energy.h>
#include <faunus/space.h>
#include <faunus/inputfile.h>
#include <faunus/geometry.h>
#include <faunus/textio.h>

namespace Faunus {

  namespace Analysis {
    AnalysisBase::AnalysisBase() : w(30), cnt(0) {
      stepcnt = 0;
    }

    AnalysisBase::AnalysisBase(Tmjson &j) : w(30), cnt(0) {
      steps = j["nstep"] | int(0);
      stepcnt = 0;
    }

    AnalysisBase::~AnalysisBase() { }

    string AnalysisBase::_info() { return string(); }

    Tmjson AnalysisBase::_json() { return Tmjson(); }

    void AnalysisBase::_sample() {
      assert(!"We should never reach here -- implement _sample() function");
      /* make pure virtual! */
    }

    void AnalysisBase::sample() {
      stepcnt++;
      if (stepcnt == steps) {
        cnt++;
        stepcnt=0;
        timer.start();
        _sample();
        timer.stop();
      }
    }

    void AnalysisBase::_test(UnitTest &t) {}

    void AnalysisBase::test(UnitTest &t) { _test(t); }

    string AnalysisBase::info() {
      using namespace textio;
      std::ostringstream o;
      if ( cnt>0 ) {
        o << header("Analysis: "+name);
        if (!cite.empty())
          o << pad(SUB,w,"Reference:") << cite << "\n";

        o << pad(SUB,w,"Sample interval") << steps << "\n";

        if (cnt>0) {
          o << pad(SUB,w,"Number of sample events") << cnt << "\n";
          double time = timer.result();
          if (time>1e-3)
            o << pad(SUB,w,"Relative time") << time << "\n";
        }
        o << _info();
      }
      return o.str();
    }

    Tmjson AnalysisBase::json() {
      Tmjson j;
      if ( ! name.empty() )
        if ( cnt>0 ) {
          j[ name ] = {
            { "samples", cnt },
            { "relative time", timer.result() }
          };
          j = merge( j, _json() );
        }
      return j;
    }

    TwobodyForce::TwobodyForce(InputMap &in, Group &g1, Group &g2, Group &_ions) {
      name="Twobody mean force calculation";
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

    void BilayerStructure::_test(UnitTest &t) {
      t("bilayer_order", S.avg() );
      t("bilayer_area", A.avg() );
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

    void CombinedAnalysis::test(UnitTest &test) {
      for (auto i : v) i->test(test);
    }

    Tmjson CombinedAnalysis::json() {
      Tmjson js;
      for (auto i : v)
        js = merge(js, i->json() );
      return js;
    }

  }//namespace
}//namespace
