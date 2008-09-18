#include "faunus/slump.h"

namespace Faunus {
  // Baseclass
  double random::random_half() { return -0.5 + random_one(); }
  bool random::runtest(float f) { return (random_one()<f) ? true : false; }

  // "slump" - the default generator
  randomDefault::randomDefault() { rand_max_inv = 1./RAND_MAX; }
  double randomDefault::random_one() { return rand_max_inv*rand(); }
  void randomDefault::random_seed(unsigned int s) { (s!=0) ? srand(s) : srand(time(0)); }

  // "Twister" - the default generator
  double randomTwister::random_one() { return mt.rand(); }
  void randomTwister::random_seed(unsigned int s) { mt.seed(s); }
}

