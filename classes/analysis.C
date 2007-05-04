#include "analysis.h"
//! \param u_init Initial total system energy
systemenergy::systemenergy(double u_init) {
  u0=u_init;
  sum=u0;
  cur=u0;
}
void systemenergy::setcurrent(double energy) { cur=energy; }
void systemenergy::operator+=(double du) { sum+=du; }
