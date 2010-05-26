#include <iostream>
#include <algorithm>
#include "faunus/average.h"
#include "faunus/xytable.h"

using namespace std;
using namespace Faunus;

int main() {
  xytable<float,int> f(0.1, -10);
  float x=-5.6;
  f(x)=10;
  cout << f(x); // --> 10

  xytable<float, average<float> > h(0.5);
  float r=7.5;
  h(r)+=2;
  h(r)+=4;
  cout << h(r).avg(); // --> 3
}
