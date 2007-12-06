#ifndef MCLOOP_H
#define MCLOOP_H

#include "inputfile.h"
#include "countdown.h"

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
    countdown<unsigned short> cnt;
  public:
    unsigned short macro, micro;
    mcloop(inputfile &);           //!< Setup
    string timing(unsigned short); //!< Show macrostep middle time and ETA.
};

/*!
 * \param in Inputfile class. "macrosteps" and "microsteps" will be searched for.
 */
mcloop::mcloop(inputfile &in)
  : cnt( in.getint("macrosteps",10))
{
  macro=in.getint("macrosteps",10);
  micro=in.getint("microsteps");
}

string mcloop::timing(unsigned short mac) {
  ostringstream o;
  o << "# Macrostep " << mac << " completed. ETA: "
    << cnt.eta(mac);
  return o.str();
}

#endif
