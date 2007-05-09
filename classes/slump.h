#ifndef _slump_h
#define _slump_h
#include <cmath>
#include <cstdlib>
#include <ctime>

/*!
 * \brief Random number functions
 * \author Mikael Lund
 * \date Lund, 2002
 */
class slump {
 private:
  double rand_max_inv;
 public:
  slump();
  void random_seed(unsigned int=0);     //!< Seed random generator (globally)
  bool runtest(float=0.5);              //!< Probability bool
  double random_one();                  //!< Random number between [0:1[
  double random_half();                 //!< Random number between [-0.5:0.5[
};
#endif

