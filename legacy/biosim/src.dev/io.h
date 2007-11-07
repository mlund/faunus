#ifndef _io_h
#define _io_h
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include "species.h"
#include "group.h"

class io {
 private:
  struct filefmt {
    double x,y,z,mw,charge;
    string atomname, aminoname;
    unsigned int atomnr, aminonr;
  };
  double charge(vector<particle> &);
  point cm(vector<particle> &);

 public:
  vector<particle> loadaam(species &, string);
  bool saveaam(species &, string, vector<particle> &, group &);
};

#endif
