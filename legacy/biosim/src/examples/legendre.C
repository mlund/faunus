#include "../legendre.h"
#include <iostream>

int main() {
  double sum=0;
  legendre l(5);

  l.eval(0.4);
  for (int i=0; i<=5; i++) {
    cout << i << " " << l.p[i] << endl;
    sum+=l.p[i];
  };
  cout << "sum = " << sum << endl;
};
