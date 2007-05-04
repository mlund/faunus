#include "analysis.h"
//! \param u_init Initial total system energy
systemenergy::systemenergy(double u_init) { u=0; u0=u_init; }
void systemenergy::operator+=(double du) { u+=du; }
