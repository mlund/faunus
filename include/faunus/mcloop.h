#ifndef FAU_MCLOOP_H
#define FAU_MCLOOP_H
#include <faunus/common.h>
#include <faunus/inputfile.h>
#include <faunus/countdown.h>

namespace Faunus {
  /*!
   * This class keeps track of the outer and inner
   * Markov chain loops. After each macro step an
   * estimated time of the simulations will be evaluated.
   *
   * The constructor will search the passed Faunus::inputfile object for the
   * keywords:
   * \li "macrosteps" - corresponding to the number of steps in the outer loop
   * \li "microsteps" - ...and the steps in the inner loop
   *
   * A typical usage is as follows\n
   * \code
   * mcloop(in) mc;
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
  class mcloop {
    private:
      countdown<unsigned int> cnt;
      bool loadstateBool;          //!< load state file if present?
      string statefile;            //!< Default name of state file to load/save
    public:
      unsigned int cnt_micro, cnt_macro;
      bool eq;
      unsigned int macro;          //!< Number of macrosteps
      unsigned int micro;          //!< Number of microsteps
      unsigned int count();        //!< Current number of steps
      mcloop(inputfile &);         //!< Setup
      string info();               //!< Get information
      string timing(unsigned int); //!< Show macrostep middle time and ETA (outdated!)
      string timing();             //!< Show macrostep middle time and ETA.
      bool macroCnt();             //!< Increase and test macro loop counter
      bool microCnt();             //!< Increase and test micro loop counter
      bool savestate(string="");   //!< Save loop state to disk
      bool loadstate(string="");   //!< Load loop state from disk
  };
}
#endif
