#include <iostream>
#include "point.h"

int main() {
  point p, o;

  p.x = 1.;
  p.y = 2.;
  p.z = 3.;

  p += p;

  cout << p << " " << o << endl;
  
};
