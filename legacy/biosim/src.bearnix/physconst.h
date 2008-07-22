#ifndef _physconst_h
#define _physconst_h

#include <iostream>

//! Placeholder for physical constants and parameters.
class physconst {
 public:
  double k,e,Na,pi,e_0,T,e_r;
  double beta,ecf,beta_ecf;
  physconst(double=298.15, double=80.);

/*
  template <class T>
    inline T ipow(T x, int n) {
      T p=1;
      for (int i=0; i<n; i++) { p*=x; };
      return p;
    };
  */
};
#endif
