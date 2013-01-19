#include <faunus/mcloop.h>
#include <faunus/inputfile.h>
#include <faunus/textio.h>

namespace Faunus {

  MCLoop::MCLoop(InputMap &in, string pfx) : cnt( in.get<int>(pfx+"macrosteps",10)) {
    string prefix=pfx;
    macro=in.get<int>(prefix+"macrosteps",10);
    micro=in.get<int>(prefix+"microsteps",0);
    cnt_micro=cnt_macro=0;
  }

  string MCLoop::info() {
    using namespace textio;
    char w=25;
    std::ostringstream o;
    o << header("MC Steps and Time")
      << pad(SUB,w,"Steps (macro micro tot)") << macro << "\u2219" << micro << " = " << macro*micro << endl
      << pad(SUB,w,"Remaining steps") << macro*micro - count() << endl;
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

}//namespace
