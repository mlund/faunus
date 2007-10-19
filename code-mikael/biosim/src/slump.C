#include "slump.h"

slump::slump() {
  rand_max_inv = 1./RAND_MAX;
};

void slump::seed(unsigned int s=13) {
  srand(s);
  //srand(time(0));
};
