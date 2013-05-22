#ifndef FAU_MCLOOP_H
#define FAU_MCLOOP_H

#ifndef SWIG
#include <faunus/common.h>
#include <faunus/textio.h>
#endif

namespace Faunus {
  class _inputfile;
  /**
   * @brief Estimate speed of a computational process
   * 
   * This class can be used to estimate when a computational
   * process will finish. The current resolution is seconds.
   */
  template<class T=unsigned int>
    class CountDown {
      private:
        T max;
        int time_i;              //!< Starting time in seconds
        std::time_t rawtime; 
        struct tm *timeinfo;
      public:
        CountDown(T);
        float speed(T);          //!< Calculate speed
        int elapsed();           //!< Time elapsed in seconds
        std::string eta(T);      //!< Estimate time of arrival
    };
  /** @param maxvalue Value at arrival */
  template<class T> CountDown<T>::CountDown(T maxvalue) {
    time_i=std::time(0);
    max = maxvalue;
  }
  template<class T> float CountDown<T>::speed(T midvalue) {
    return float(std::time(0)-time_i) / midvalue;
  }
  /** 
   * @param midvalue Value somewhere between start and arrival
   * @return String with estimated time and date of arrival
   */
  template<class T> std::string CountDown<T>::eta(T midvalue) {
    rawtime = std::time(NULL) + int( speed(midvalue) * (max-midvalue)  );
    timeinfo = std::localtime(&rawtime);
    return std::asctime(timeinfo);
  }
  template<class T> int CountDown<T>::elapsed() { return time(0)-time_i; }

  /**
   * @brief Two level loop book-keeping
   *
   * This class simply keeps track of an outer and inner
   * Markov chain loop. After each macro step an
   * estimated time of the simulation will be evaluated.
   *
   * The constructor will search the passed Faunus::InputMap object for the
   * keywords:
   *
   * Key               | Description
   * :---------------- | :-----------------------------
   * `loop_macrosteps` | Number of steps in outer loop
   * `loop_microsteps` | Number of steps in inner loop
   *
   * Typical usage:
   *
   * @code
   *
   * InputMap mcp("myinput");
   * MCLoop(mcp) mc;
   *
   * while ( mc.macroCnt() ) {
   *   while ( mc.microCnt() ) {
   *   }
   *   std::cout << mc.timing();
   * }
   * std::cout << mc.info();
   *
   * @endcode
   * @date 2007
   * @todo This could be made general with an arbitrary number of levels
   */
  class MCLoop {
    private:
      CountDown<unsigned int> cnt;
      unsigned int macro;          //!< Number of macrosteps
      unsigned int micro;          //!< Number of microsteps
      unsigned int cnt_micro, cnt_macro;
      bool eq;
      string timing(unsigned int); //!< Show macrostep middle time and ETA (outdated!)
    public:
      MCLoop(InputMap&, string="loop_"); //!< Setup
      unsigned int count();        //!< Current number of steps
      string info();               //!< Get information
      string timing();             //!< Show macrostep middle time and ETA.
      bool macroCnt();             //!< Increase and test macro loop counter
      bool microCnt();             //!< Increase and test micro loop counter
  };

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

}
#endif
