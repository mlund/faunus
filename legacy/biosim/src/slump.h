#ifndef _slump_h
#define _slump_h

// SLUMP CLASS - generates random numbers
// M. Lund, 2002
#include <cmath>
#include <cstdlib>
#include <ctime>

using namespace std;

class slump {
 private:
  double rand_max_inv;
 public:
  slump();
  void seed(unsigned int);       ///! Seed random generator (globally)
  double random_one() {   ///! Random number between [0:1[
    return rand_max_inv*rand();
  };

  double random_half() {  ///! Random number between [-0.5:0.5[
    return -0.5 + rand_max_inv * rand();
  };

  // Violate C++ standard!
  // Work-around for GCC: -fno-strict-aliasing
  inline float InvSqrtE(float x) {
    float xhalf = 0.5f * x;
    int i = *(int*)&x; // store floating-point bits in integer
    i = 0x5f3759d5 - (i >> 1); // initial guess for Newton's method
    x = *(float*)&i; // convert new bits into float
    x = x*(1.5f - xhalf*x*x); // 1st round of Newton's method
    x = x*(1.5f - xhalf*x*x); // 2nd round of Newton's method
    return x;
  };

};

#endif
