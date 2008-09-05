#ifndef FAU_HARDSPHERE_H
#define FAU_HARDSPHERE_H

#include <vector>
#include "point.h"
#include "group.h"

namespace Faunus {

/*!\brief Hardsphere overlap between particles
 * \author Mikael Lund
 */
class hardsphere {
  public:
    bool overlap(vector<particle> &, int);                ///< all<->particle i.
    bool overlap(vector<particle> &, particle &);         ///< all<->arbitrary (external) particle.
    bool overlap(vector<particle> &, group &);            ///< all<->group.
    bool overlap(vector<particle> &, group &, int);       ///< group<->particle i.
    bool overlap(vector<particle> &, group &, particle &);///< group<>arbitrary (external) particle
    bool overlap(vector<particle> &, group &, group &);   ///< group<->group.
    bool overlap(vector<particle> &, group::group &, group::group &, double &);   ///< group<->group.
    bool overlap(vector<particle> &, vector<short int> &, double);      ///< internal collisions within subset
};
}
#endif
