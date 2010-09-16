#ifndef FAU_REPLICAEXCHANGE_H
#define FAU_REPLICAEXCHANGE_H

#include "faunus/common.h"
#include "faunus/vectorid.h"
#include "faunus/ensemble.h"

namespace Faunus {

  class simBottle;

  /*!
   * \brief Move for generic replica exchange
   * \author Mikael Lund
   * \date Malmo, 2010
   * \todo Optimize by swapping only particle coordinates.
   *
   * This class will try to swap configurations between replica
   * simulations (simBottle base class). Old energies are
   * taken from the current energy of the replicas, so make sure
   * they are always in sync. New energies are calculated by calling
   * the systemEnergy() function in simBottle.
   */
  class replicaexchange {
    private:
      canonical nvt;
      vectorid<int> cnt, accepted;
    public:
      bool swap(simBottle &, simBottle &);
      string info();
  };

}//namespace
#endif

