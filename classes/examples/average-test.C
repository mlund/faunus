#include <iostream>
#include "../average.h"

using namespace std;

int main() {
  average<float> x;
  x+=2.0;
  x+=6.0;
  cout << x.avg(); // --> 4.0
}
