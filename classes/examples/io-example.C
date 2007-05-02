#include <iostream>
#include <iterator>
#include "io.h"

using namespace std;

int main() {
  species spc;
  ioaam aam(spc);
  vector<particle> p=aam.load("io-example.dat");
  aam.save("io-example.out", p);

  iopov pov(spc);
  pov.save("io-example.pov", p); 
  
}
