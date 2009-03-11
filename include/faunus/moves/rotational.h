#ifndef FAU_ROTATIONAL_H
#define FAU_ROTATIONAL_H

#include "faunus/moves/base.h"

namespace Faunus {
  
  /*! \brief Rotate group around its mass-center.
   *  \author Mikael Lund
   *  \date Prague 2007
   */
  class macrorot : public markovmove { 
    public:
      macrorot( ensemble&, container&, energybase&);
      double move(macromolecule &);
  };

}

#endif
