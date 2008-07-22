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

  inline double random_one() {   ///! Random number between [0:1[
    return rand_max_inv*rand();
  };

  inline double random_half() {  ///! Random number between [-0.5:0.5[
    return -0.5 + rand_max_inv * rand();
  };
};
#endif
