#ifndef _chain_h
#define _chain_h

#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include "species.h"
#include "point.h"

using namespace std;

class chain {
 public:
  int graftpoint;        //Graftpoint in particle vector. -1 if not grafted.
  vector<particle> v;
  chain(species &,string, int=-1); //residue string , graftpoint 
};

#endif
