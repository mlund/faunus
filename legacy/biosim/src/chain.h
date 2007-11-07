#ifndef _chain_h
#define _chain_h

#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include "peptide.h"
#include "point.h"

using namespace std;

class chain : private aminoacid {
 public:
  int graftpoint;        //Graftpoint in particle vector. -1 if not grafted.
  vector<particle> v;
  chain(string, int=-1); //residue string , graftpoint 
};

#endif
