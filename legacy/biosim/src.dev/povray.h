#ifndef _povray_h
#define _povray_h

#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include "point.h"
#include "group.h"

using namespace std;

class povray {
 public:
  string anionTxt, cationTxt, neutralTxt;
  ostringstream ostr;

  povray();
  void zaxis(double);     //create z-axis from [-double:double]
  void bond(point &, point &, float=1.);
  void coordinates(float); //print axes of coord. system
  void sphere(double); //spherical environment
  void glowingsphere(particle &);
  void cube(double); //quadratic box
  void cylinder(point &, point &, double=0.5);
  void add(vector<particle> &, group &);
  void save(string); //give filename

 private:
  void defineTextures();
};

#endif
