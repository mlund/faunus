#ifndef FAU_MCLOOP_H
#define FAU_MCLOOP_H

#ifndef SWIG
#include <faunus/common.h>
#include <faunus/textio.h>
#endif

namespace Faunus {
  class _inputfile;
  /*
   * \brief Estimate speed of a computational process
   * \author Mikael Lund
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
   * \brief Monte Carlo loop book-keeping
   *
   * This class keeps track of the outer and inner
   * Markov chain loops. After each macro step an
   * estimated time of the simulations will be evaluated.
   *
   * The constructor will search the passed Faunus::InputMap object for the
   * keywords:
   * \li "macrosteps" - corresponding to the number of steps in the outer loop
   * \li "microsteps" - ...and the steps in the inner loop
   *
   * A typical usage is as follows\n
   * \code
   * MCLoop(in) mc;
   * while (mc.macroCnt()) {
   *   while (mc.microCnt()) {
   *     //inner loop code
   *   }
   *   std::cout << mc.timing();
   * }
   * std::cout << mc.info();
   * \endcode
   *
   * \author Mikael Lund
   * \date 2007
   */
  class MCLoop {
    private:
      CountDown<unsigned int> cnt;
      bool loadstateBool;          //!< load state file if present?
      string statefile;            //!< Default name of state file to load/save
      bool savestate(string="");   //!< Save loop state to disk
      bool loadstate(string="");   //!< Load loop state from disk
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
