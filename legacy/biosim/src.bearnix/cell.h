#ifndef _cell_h
#define _cell_h

#include <vector>
#include "slump.h"
#include "point.h"
#include "physconst.h"

//using namespace std;

class cell : private physconst, private slump {
 private:
  double vol_A, vol_L;

 public:
  cell(double r);
  double cell_r, cell_r2, cell_dia;
  inline void randomPos( point &p );
  inline bool cellCollision( particle &p );
  double conc(int);
  double ionStr( vector<particle> & );
  inline double segmentVolume(double,double);
};

// generate random position within cell
// (this could be made more elegant using polar coords)
inline void cell::randomPos( point &p ) {
  double l=cell_r2+1;
  while (l>cell_r2) {
    p.x = random_half() * cell_dia ;
    p.y = random_half() * cell_dia ;
    p.z = random_half() * cell_dia ;
    l=p.x*p.x+p.y*p.y+p.z*p.z; //squared distance from origo
  };
};
// collision check with cell wall
inline bool cell::cellCollision( particle &p ) {
  if (p.x*p.x+p.y*p.y+p.z*p.z > cell_r2) return true;
  return false;
};


#endif
