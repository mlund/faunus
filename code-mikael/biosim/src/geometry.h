/*
 * GEOMETRY (M.Lund 2002)
 *
 * Calculate distances between
 * coordinates
 *
 */

#ifndef geometry_h
#define geometry_h

#include <cmath>
#include "point.h"
#include "intrinsics.h"

using namespace std;

class geometry {
  public:
    inline double dist(point &p1, point &p2) {
      return sqrt(sqdist(p1,p2));
    };
    inline double sqdist(point::point &p1, point::point &p2) {
      double dx,dy,dz;
      dx=p1.x-p2.x;
      dy=p1.y-p2.y;
      dz=p1.z-p2.z;
      return dx*dx + dy*dy + dz*dz;
    };
    inline double invdist(point::point &p1, point::point &p2) {
      return 1./sqrt( sqdist(p1, p2) );
    };
};
#endif

