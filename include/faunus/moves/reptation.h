#ifndef FAU_REPTATION_H
#define FAU_REPTATION_H

#include "faunus/moves/base.h"
#include "faunus/energy/springinteraction.h"

namespace Faunus {

  /*!
   * Reptation move for polymers
   *
   * \author Bjoern Persson
   * \date Lund, 2010
   */
  class reptation : public markovmove {
    private:
      enum keywords {FORWARD, BACKWARD};
      keywords direction;
      unsigned short beg;            //!< First CS fix point
      unsigned short end;            //!< Last CS fix point
      unsigned short len;            //!< Number of atoms to include in crank-shaft move

    public:

      reptation(ensemble &, container &, energybase &, inputfile &);
      double move(polymer &);        //!< Move, assuming rigid, equidistant bonds.
      string info();
  };
} // namespace
#endif
