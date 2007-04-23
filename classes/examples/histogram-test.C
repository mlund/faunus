#include <iostream>
#include "../histogram.h"

using namespace std;

int main() {
  histogram h(0.5);
  h.add(0.5);
  h.add(1.0);
  h.add(1.0);
  h.show();
  float val=h.get(0.5); // val=0.3333
}

// Output:
// 0.5 0.333333
// 1 0.666667
