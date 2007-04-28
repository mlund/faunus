#ifndef _metropolis_h
#define _metropolis_h
#include "slump.h"

//! Polymorph ensemble class
class ensemble {
  protected:
    slump slp;
  public:
    virtual bool metropolis(double)=0; //!< Accept critetia
};

/*! \brief NVT ensemble
 */
class canonical : public ensemble {
  public:
    bool metropolis(double du)
    {
      if (du > 0)
        if ( slp.random_one()>exp(-du) )
          return false;
      return true;
    }
};

class grandcanonical : public ensemble {};
#endif
