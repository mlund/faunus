#include <iostream>
#include "../average.h"
#include "../xytable.h"

using namespace std;

int main() {
  xytable<float,int> f(0.1, -10);
  float x=-5.6;
  f(x)=10;
  cout << f(x) << endl;

  xytable<float, average<float> > h(0.5);
  float r=7.5;
  h(r)+=2;
  h(r)+=4;
  cout << h(r).avg() << endl; 
}
