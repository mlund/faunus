/*  CELL CLASS
    -------------------------------------------*/
#include "cell.h"

cell::cell(double r) {
  cell_r = r;          //cell radius
  cell_r2 = r*r;       //squared radius
  cell_dia = 2*r;      //diameter
  vol_A = (4./3.)*pi*r*r*r;//volume in AA^3
  vol_L = vol_A*1e-27; //volume in liters
};

//return molar concentration
double cell::conc(int number) {
  return (number/Na) / vol_L;
};


//calculate volume of spherical segment, starting
//d from the centrum, and of thickness h.
//- see Wolfram research: "Spherical Segment"
double cell::segmentVolume(double d, double h) {
  return pi*h*( cell_r2 - d*d - h*d - h*h/3. );
};


double cell::ionStr(vector<particle> &p) {
  return 0;
};
