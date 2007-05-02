/*  SIMBOX CLASS
    -------------------------------------------*/
#include "simbox.h"

Simbox::Simbox(double boxlen, double boxsep) {
  len = boxlen;             //lateral boxlength
  sep = boxsep;             //sepration between surfaces
  area = len*len;        //area of a single surfaces
  len_half = len*0.5;    //half length
  sep_half = sep*0.5;    //half separation
  sqlen_half = sep_half*sep_half;
  vol_A = len*len*sep;          //volume in AA^3
  vol_L = vol_A*1e-27;          //volume in liters
};

//return molar concentration
double Simbox::conc(int number) {
  return (number/Na) / vol_L;
};

//calculate volume of spherical segment, starting
//d from the centrum, and of thickness h.
//- see Wolfram research: "Spherical Segment"
//double simbox::segmentVolume(double d, double h) {
//  return pi*h*( cell_r2 - d*d - h*d - h*h/3. );
//};


//double simbox::ionStr(vector<particle> &p) {
//  return 0;
//};
