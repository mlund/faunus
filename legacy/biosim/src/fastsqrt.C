#include <iostream>
#include <cmath>
#include "slump.h"

using namespace std;

int main() {
  slump f;
  for (float x=0.1; x<200; x+=0.1) {
    double a=1./f.InvSqrtE(x*x);
    double b=sqrt(x*x);
    cout << a - b << endl;
  };
};
