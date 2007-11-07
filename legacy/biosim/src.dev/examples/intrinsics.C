#include <iostream>
#include <cmath>
#include "../intrinsics.h"

using namespace std;

int main() {

  cout << "frsqrte:" << endl
       << "  manual = " << 1./sqrt(25.) << endl
       << "  intrinsic = " << frsqrte(25.) << endl;
  
};
