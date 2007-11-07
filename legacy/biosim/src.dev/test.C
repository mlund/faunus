#include "hist.h"
#include <limits>
#include <iostream>

using namespace std;

int main() {
  distribution<double> h(0.5, "g(r)");
  h.add(1,2);
  h.add(2,3);
  h.add(2,4);
  cout << h.show();
};
