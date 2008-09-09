#ifndef FAU_METROPOLIS_H
#define FAU_METROPOLIS_H

#include "faunus/slump.h"

namespace Faunus {
/*!
 * Polymorph ensemble class
 * \todo Maybe not really needed...
 */
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

}//namespace
#endif
