#include "slump.h"
namespace Faunus {
  slump::slump() { rand_max_inv = 1./RAND_MAX; }
  double slump::random_one() { return rand_max_inv*rand(); }
  double slump::random_half() { return -0.5 + rand_max_inv * rand(); }
  void slump::random_seed(unsigned int s) { (s!=0) ? srand(s) : srand(time(0)); }
  bool slump::runtest(float f) { return (random_one()<f) ? true:false; }
}
