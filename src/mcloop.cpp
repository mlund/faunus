#include <faunus/mcloop.h>

namespace Faunus {

  mcloop::mcloop(inputfile &in) : cnt( in.getint("macrosteps",10)) {
    macro=in.getint("macrosteps",10);
    micro=in.getint("microsteps");
    //eq=in.getboo("equilibration", false);
    cnt_micro=cnt_macro=0;
  }

  string mcloop::info() {
    std::ostringstream o;
    o << endl << "# STEP AND TIME DETAILS:" << endl
      << "#   Steps (macro micro tot)= "
      << macro << " x " << micro << " = " << macro*micro << endl;
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
   * maximum value has been reached.
   */
  bool mcloop::macroCnt() {
    return (++cnt_macro>macro) ? false : true;
  }

  /*!
   * As mcloop::macroCnt() but for the microsteps
   */
  bool mcloop::microCnt() {
    if (cnt_micro++<micro)
      return true;
    cnt_micro=0;
    return false;
  }

}//namespace
