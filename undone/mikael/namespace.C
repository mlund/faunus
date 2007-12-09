#include "point.h"
#include <iostream>

namespace faunus {
  double abe;
  class test(double) { return; };
}

namespace faunus {
  double kat;
  faunus::test t(2.0);
}

using namespace faunus;
using namespace std;

int main() {
  abe=0.1;
  kat=0.2;
  point a;
  cout << abe << " " << kat << endl;
  return 0;
};
