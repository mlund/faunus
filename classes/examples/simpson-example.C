#include <iostream>
#include "../simpson.h"

using namespace std;

double function(double x) { return x*x; }

int main() {
  cout << simpson<function>(0,3,100) << endl;
};
