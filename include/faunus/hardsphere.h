#ifndef FAU_HARDSPHERE_H
#define FAU_HARDSPHERE_H
#include "faunus/point.h"
#include "faunus/group.h"

namespace Faunus {
  /*!\brief Hardsphere overlap between particles
   * \author mikaek lund
   */
  class hardsphere {
    public:
      bool overlap(const vector<particle> &, int);                            ///< all<->particle i.
      bool overlap(const vector<particle> &, const particle &);               ///< all<->arbitrary (external) particle.
      bool overlap(const vector<particle> &, const group &);                  ///< all<->group.
      bool overlap(const vector<particle> &, const group &, int);             ///< group<->particle i.
      bool overlap(const vector<particle> &, const group &, const particle &);///< group<>arbitrary (external) particle
      bool overlap(const vector<particle> &, const group &, const group &);   ///< group<->group.
      bool overlap(const vector<particle> &, const group &, const group &, double &); ///< group<->group.
      //bool overlap(vector<particle> &, vector<short int> &, double);          ///< internal collisions within subset
  };
}
#endif
