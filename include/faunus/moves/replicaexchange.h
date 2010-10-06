#ifndef FAU_REPLICAEXCHANGE_H
#define FAU_REPLICAEXCHANGE_H

#include "faunus/common.h"
#include "faunus/vectorid.h"
#include "faunus/ensemble.h"
#include "faunus/average.h"

namespace Faunus {

  class bottle;

  /*!
   * \brief Move for generic replica exchange
   * \author Mikael Lund
   * \date Malmo, 2010
   * \todo Optimize by swapping only particle coordinates.
   *
   * This class will try to swap configurations between replica
   * simulations (bottle base class). Old energies are
   * taken from the current energy of the replicas, so make sure
   * they are always in sync. New energies are calculated by calling
   * the bottle::systemEnergy() function.
   */
  class replicaexchange {
    private:
      canonical nvt;
      vectorid<int> cnt, accepted;
      vectorid< average<double> > dU, expdU;
    public:
      bool swap(bottle &, bottle &);
      string info();
  };

}//namespace
#endif

