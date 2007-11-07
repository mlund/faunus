#ifndef _physconst_h
#define _physconst_h

#include <iostream>

//! Placeholder for physical constants and parameters.
class physconst {
 public:
  double k,e,Na,pi,e_0,T,e_r;
  double beta,ecf,beta_ecf;
  physconst(double=298.15, double=80.);
};
#endif
