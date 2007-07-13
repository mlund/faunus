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
