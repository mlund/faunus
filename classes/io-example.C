#include <iostream>
#include <iterator>
#include "io.h"

using namespace std;

int main() {
  vector<string> v;
  io i;
  i.readfile("io-example.dat", v);
  i.strip(v);
  for (int i=0; i<v.size(); i++)
    cout << v[i] << endl;

  species spc;
  ioaam aam(spc);
  vector<particle> p=aam.load("io-example.dat");
  for (int i=0; i<p.size(); i++)
    cout << p[i] << endl;
  aam.save(p, "io-example.out");

  //iopov pov(spc);
  //pov.save(p, "io-example.pov"); 
  
}
