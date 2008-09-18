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
   * \author mikaek lund
   * \date 2007
   */
  class mcloop {
    private:
      countdown<unsigned int> cnt;
    public:
      unsigned int macro, micro;
      bool eq;
      mcloop(inputfile &);         //!< Setup
      string info();               //!< Shows setup
      string timing(unsigned int); //!< Show macrostep middle time and ETA.
  };
  /*!
   * \param in Inputfile class. "macrosteps" and "microsteps" will be searched for.
   */
  mcloop::mcloop(inputfile &in)
    : cnt( in.getint("macrosteps",10))
  {
    macro=in.getint("macrosteps",10);
    micro=in.getint("microsteps");
    eq=in.getboo("equilibration", false);
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
  string mcloop::timing(unsigned int mac) {
    std::ostringstream o;
    o << "# Macrostep " << mac << " completed. ETA: "
      << cnt.eta(mac);
    return o.str();
  }
}
#endif
