#ifndef _LENNARDJONES_H
#define _LENNARDJONES_H

#include "point.h"

/*!
 *  \brief Lennard-Jones potential
 *  \author Mikael Lund
 *  \date Prague, 2007
 */
class pot_lj {
  public:
    double eps;
    pot_lj(double epsilon) { eps=epsilon; }
    /*!
     *  L-J pair energy.
     *  \f$ u_{lj} = \epsilon \left ( \frac{\sigma}{r^{12}} - \frac{\sigma}{r^6} \right ) \f$
     *  \param r2 Squared distance between particle 1 and 2.
     */
    inline double lj(particle &p1, particle &p2, double &r2) {
      register double x=p1.radius+p2.radius,
             u=x*x/r2;
      x=u*u*u;
      return (x*x-x)*eps;
    }
};

#endif
