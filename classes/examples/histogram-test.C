#include <iostream>
#include "../histogram.h"

using namespace std;

int main() {
  histogram h(0.5);
  h.add(0.5);
  h.add(1.0);
  h.add(1.0);

  cout << h.get(0.5) << " "
       << h.get(1.0) << endl;
}
