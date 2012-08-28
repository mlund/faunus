#include <faunus/mcloop.h>
#include <faunus/inputfile.h>
#include <faunus/textio.h>

namespace Faunus {

  MCLoop::MCLoop(InputMap &in, string pfx) : cnt( in.get<int>(pfx+"macrosteps",10)) {
    string prefix=pfx;
    macro=in.get<int>(prefix+"macrosteps",10);
    micro=in.get<int>(prefix+"microsteps",0);
    statefile=in.get<string>(prefix+"statefile", "loop.state");
    loadstateBool = false; //in.get<bool>(prefix+"loadstate", false);
    if (loadstateBool)
      loadstate();
    //eq=in.get<bool>("equilibration", false);
    cnt_micro=cnt_macro=0;
  }

  string MCLoop::info() {
    using namespace textio;
    char w=25;
    std::ostringstream o;
    o << header("MC Steps and Time")
      << pad(SUB,w,"Steps (macro micro tot)") << macro << "\u2219" << micro << " = " << macro*micro << endl
      << pad(SUB,w,"Remaining steps") << macro*micro - count() << endl;
    if (loadstateBool)
      o << "#   Load state from disk   = yes" << endl;
    int t=cnt.elapsed();
    if (t>5) {
      o << pad(SUB,w,"Time elapsed") << t/(3600.) << " h" << endl
        << pad(SUB,w,"Steps/minute") << macro*micro/(t/60.) << endl;
    }
    o << std::flush;
    return o.str();
  }

  string MCLoop::timing() {
    return timing(cnt_macro);
  }

  /*!
   * \note This will try to flush the output stream buffer.
   */
  string MCLoop::timing(unsigned int mac) {
    using namespace textio;
    std::ostringstream o;
    o << indent(SUB) << "Macrostep " << std::left << std::setw(5) << mac << "ETA: "
      << cnt.eta(mac) << std::flush;
    return o.str();
  }

  /*!
   * Increase macroloop counter and test if the
   * maximum value has been reached. Whenever called
   * this function saves a state file to disk
   */
  bool MCLoop::macroCnt() {
    //savestate();
    return (++cnt_macro>macro) ? false : true;
  }

  /*!
   * As MCLoop::macroCnt() but for the microsteps
   */
  bool MCLoop::microCnt() {
    if (cnt_micro++<micro)
      return true;
    cnt_micro=0;
    return false;
  }

  /*!
   * Returns the number of completed steps
   */
  unsigned int MCLoop::count() {
    return (cnt_macro-1)*micro + cnt_micro;
  }

  bool MCLoop::savestate(string name) {
    if (name.empty())
      name=statefile;
    std::ofstream f(name.c_str());
    if (f) {
      f << macro << " " << micro << " " << cnt_macro << " " << cnt_micro << " " << count();
      f.close();
      return true;
    }
    return false;
  }

  bool MCLoop::loadstate(string name) {
    unsigned int _macro, _micro, _cnt_macro, _cnt_micro, _cnt;
    if (loadstateBool) {
      if (name.empty())
        name=statefile;
      std::ifstream f(name.c_str());
      if (f) {
        f >> _macro >> _micro >> _cnt_macro >> _cnt_micro >> _cnt;
        f.close();
        // Simple extension if micro/macro unchanged
        if (_micro==micro && _macro==macro) {
          cnt_micro=_cnt_micro;
          cnt_macro=_cnt_macro;
          macro+=cnt_macro;
          return true;
        }
        // Otherwise map onto current micro/macro values
        cnt_macro = _cnt / micro;
        cnt_micro = _cnt - (cnt_macro)*micro;
        macro+=cnt_macro;
        return true;
      }
    }
    return false;
  }
}//namespace
