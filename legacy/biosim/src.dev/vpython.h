#ifndef _vpython_h
#define _vpython_h

#include "point.h"
#include "group.h"
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>

class vpython {
  private:
    ostringstream ostr;
  public:
    vpython();
    void add(vector<particle> &, group &);
    bool save(string);
};

class vrml {
  private:
    ostringstream ostr;
  public:
    vrml();
    void add(vector<particle> &, group &);
    bool save(string);
};

#endif
