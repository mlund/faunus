#ifndef FAU_MCLOOP_H
#define FAU_MCLOOP_H

#ifndef SWIG
#include <faunus/common.h>
#include <faunus/textio.h>
#endif

namespace Faunus {
  class _inputfile;
  /*
   * @brief Estimate speed of a computational process
   * 
   * This class can be used to estimate when a computational
   * process will finish. The current resolution is seconds.
   */
  template<class T>
    class CountDown {
      private:
        T max;
        int time_i;              //!< Starting time in seconds
        time_t rawtime; 
        struct tm *timeinfo;
      public:
        CountDown(T);
        float speed(T);          //!< Calculate speed
        int elapsed();           //!< Time elapsed in seconds
        std::string eta(T);      //!< Estimate time of arrival
    };
  /*! \param maxvalue Value at arrival
  */
  template<class T> CountDown<T>::CountDown(T maxvalue) {
    time_i=time(0);
    max=maxvalue;
  }
  template<class T> float CountDown<T>::speed(T midvalue) {
    return float(time(0)-time_i) / midvalue;
  }
  /*! \param midvalue Value somewhere between start and arrival
   *  \return String with estimated time and date of arrival
   */
  template<class T> std::string CountDown<T>::eta(T midvalue) {
    rawtime = time(NULL) + int( speed(midvalue) * (max-midvalue)  );
    timeinfo = localtime(&rawtime);
    return asctime(timeinfo);
  }
  template<class T> int CountDown<T>::elapsed() { return time(0)-time_i; }

  /*!
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
}
#endif
