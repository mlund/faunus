#include <faunus/mcloop.h>

namespace Faunus {

  mcloop::mcloop(inputfile &in) : cnt( in.getint("macrosteps",10)) {
    string prefix="loop_";
    macro=in.getint("macrosteps",10);
    micro=in.getint("microsteps");
    statefile=in.getstr(prefix+"statefile", "loop.state");
    loadstateBool = false; //in.getboo(prefix+"loadstate", false);
    if (loadstateBool)
      loadstate();
    //eq=in.getboo("equilibration", false);
    cnt_micro=cnt_macro=0;
  }

  string mcloop::info() {
    std::ostringstream o;
    o << endl << "# STEP AND TIME DETAILS:" << endl
      << "#   Steps (macro micro tot)= " << macro << " x " << micro << " = " << macro*micro << endl
      << "#   Remaining steps        = " << macro*micro - count() << endl;
    if (loadstateBool)
      o << "#   Load state from disk   = yes" << endl;
    int t=cnt.elapsed();
    if (t>5) {
      o << "#   Time elapsed (hours)   = " << t/(3600.) << endl
        << "#   Steps/minute           = " << macro*micro/(t/60.) << endl;
    }
    return o.str();
  }

  string mcloop::timing() {
    return timing(cnt_macro);
  }

  /*!
   * \note This will try to flush the output stream buffer.
   */
  string mcloop::timing(unsigned int mac) {
    std::ostringstream o;
    o << "# Macrostep " << mac << " completed. ETA: "
      << cnt.eta(mac) << std::flush;
    return o.str();
  }

  /*!
   * Increase macroloop counter and test if the
   * maximum value has been reached. Whenever called
   * this function saves a state file to disk
   */
  bool mcloop::macroCnt() {
    savestate();
    return (++cnt_macro>macro) ? false : true;
  }

  /*!
   * As mcloop::macroCnt() but for the microsteps
   */
  bool mcloop::microCnt() {
    savestate();
    if (cnt_micro++<micro)
      return true;
    cnt_micro=0;
    return false;
  }

  /*!
   * Returns the number of completed steps
   */
  unsigned int mcloop::count() {
    return (cnt_macro-1)*micro + cnt_micro;
  }

  bool mcloop::savestate(string name) {
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

  bool mcloop::loadstate(string name) {
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
