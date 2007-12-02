#include "analysis.h"

//! \param u_init Initial total system energy
systemenergy::systemenergy(double u_init) {
  u0=u_init;
  sum=u0;
  cur=u0;
}
void systemenergy::update(double energy) {
  cur=energy;
  uavg+=cur;
  u2avg+=cur*cur;
}
void systemenergy::operator+=(double du) { sum+=du; }


distributions::distributions(unsigned short n, float min, float max) {
  // Resize vectors!
  // (cumbersome syntax, but there's no way around it)
  d.resize( n, xytable<float, average<float> >(0.5,min,max));
  s.resize( n );
}

void distributions::add( unsigned short i, float x, float y ) { d[i](x)+=y; }

string distributions::info() {
  
}
