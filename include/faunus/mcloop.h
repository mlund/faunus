#ifndef FAU_MCLOOP_H
#define FAU_MCLOOP_H
#include "faunus/inputfile.h"
#include "faunus/countdown.h"

namespace Faunus {
  /*!
   * This class keeps track of the outer and inner
   * Markov chain loops. Also, after each macro step
   * it can estimate when the simulation will finish
   *
   * \author Mikael Lund
   * \date 2007
   */
  class mcloop {
    private:
      countdown<unsigned int> cnt;
      unsigned int cnt_micro, cnt_macro;
    public:
      unsigned int macro, micro;
      bool eq;
      mcloop(const inputfile &);   //!< Setup
      string info();               //!< Shows setup
      string timing(unsigned int); //!< Show macrostep middle time and ETA.
      string timing();             //!< Show macrostep middle time and ETA.
      bool macroCnt();             //!< Increase macro loop counter
      bool microCnt();             //!< Increase micro loop counter
  };
  /*!
   * \param in Inputfile class. "macrosteps" and "microsteps" will be searched for.
   */
  mcloop::mcloop(const inputfile &in) : cnt( in.getint("macrosteps",10)) {
    macro=in.getint("macrosteps",10);
    micro=in.getint("microsteps");
    eq=in.getboo("equilibration", false);
    cnt_micro=cnt_macro=0;
  }
  string mcloop::info() {
    std::ostringstream o;
    o << endl << "# STEP AND TIME DETAILS:" << endl
      << "#   Steps (macro micro tot)= "
      << macro << " " << micro << " " << macro*micro << endl;
    int t=cnt.elapsed();
    if (t>5) {
      o << "#   Time elapsed (hours)   = " << t/(3600.) << endl
        << "#   Steps/minute           = " << macro*micro/(t/60.) << endl;
    }
    return o.str();
  }
  string mcloop::timing() { return timing(cnt_macro); }
  string mcloop::timing(unsigned int mac) {
    std::ostringstream o;
    o << "# Macrostep " << mac << " completed. ETA: "
      << cnt.eta(mac);
    return o.str();
  }

  bool mcloop::macroCnt() { return (++cnt_macro>macro) ? false : true; }
  bool mcloop::microCnt() {
    if (cnt_micro++<micro)
      return true;
    cnt_micro=0;
    return false;
  }
}
#endif
