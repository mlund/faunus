#include "container.h"

namespace Faunus {
//----------- CONTAINER -----------------
string container::info() {
  float z=charge();
  ostringstream o;
  o << endl
    << "# SIMULATION CONTAINER:" << endl
    << "#   Number of particles  = " << p.size() << endl
    << "#   Volume (AA^3)        = " << volume << endl
    << "#   Electroneutrality    = " 
    << ((z!=0) ? "NO!" : "Yes") << " "  << z << endl;;
  return o.str();
}
string container::povray() {
  return string(
        "#declare cell=texture {\n"
        "pigment {color rgbf <1,1,1,.9>}\n"
        " finish {phong .1 ambient .2}\n"
        "}\n" );
}

//----------- CELL ----------------------
cell::cell(float radius) { setradius(radius); }
cell::cell(inputfile &in) { setradius(in.getflt("cellradius"));}
void cell::setradius(float radius) {
  r = radius; 
  r2 = r*r; 
  diameter = 2*r; 
  volume = (4./3.)*acos(-1.)*r*r*r;
}
string cell::info() {
  ostringstream o;
  o << container::info() 
    << "#   Shape                = Spherical" << endl
    << "#   Radius               = " << r << endl;
  return o.str();
}
void cell::randompos(point &p) {
  double l=r2+1;
  while (l>r2) {
    p.x = slp.random_half()*diameter;
    p.y = slp.random_half()*diameter;
    p.z = slp.random_half()*diameter;
    l=p.x*p.x+p.y*p.y+p.z*p.z;
  }
}
string cell::povray() {
  ostringstream o;
  o << "sphere {<0,0,0>," << r << " texture {cell}}\n"
    << "cylinder {<0,0," << -r << ">,<0,0," << r << ">,0.5 texture {cell}}\n";
  return o.str();
}
//----------- BOX --------------------------
void box::setvolume(double v) {
  len=pow(v, 1./3);
  setlen(len);
}
void box::setlen(double l) {
  len = l;    // cubic box sidelength
  len_half=len/2;
  len_inv=1./len;
  volume = len*len*len;
}
box::box(double l) { setlen(l); }
box::box(inputfile &in) { setlen( in.getflt("boxlen") ); }
string box::info() {
  ostringstream o;
  o << container::info() 
    << "#   Shape                = Cube" << endl
    << "#   Side length          = " << len << endl;
  return o.str();
}
point box::randompos() {
  point p;
  randompos(p);
  return p;
}
void box::randompos(point &p) {
  p.x = slp.random_half()*len;
  p.y = slp.random_half()*len;
  p.z = slp.random_half()*len;
}
/*void box::randompos(vector<point> &p) {
  short unsigned size=p.size;
}*/
string box::povray() {
  ostringstream o;
  o << "box {<" <<-len_half <<"," <<-len_half <<"," <<-len_half <<"> , <"
    << len_half <<"," <<len_half <<"," <<len_half <<"> texture {cell}}\n";
  return o.str();
}
//----------- SLIT --------------------------
string slit::info() {
  ostringstream o;
  o << container::info() 
    << "#   Shape                = Cube - xy periodicity, only" << endl
    << "#   Side length          = " << len << endl;
  return o.str();
}

//-------------- CLUTCH -------------------------
//! \param radius Radius of the cell
//! \param min Beginning of the particle-excluded region
//! \param max Ending of the particle-excluded region (Note max>min)
clutch::clutch(float radius, float min, float max) {
  r=radius;
  r2=r*r;
  diameter=2*r;
  volume=(4./3.)*acos(-1.)*r2*r;
  zmin=min;
  zmax=max;
}
void clutch::randompos(point &p) {
  double l=r2+1;
  while (l>r2) {
    p.x = slp.random_half()*diameter;
    p.y = slp.random_half()*diameter;
    p.z = slp.random_half()*diameter;
    if (p.z>zmax || p.z<zmin)
      l=p.x*p.x+p.y*p.y+p.z*p.z; //squared distance from origo
  };
}
//------------CYLINDER---------------------------
//! \param len  Length of the cylinder
//! \param r    Radius of the cylinder
cylinder::cylinder(float length, float radius) {
  len=length;
  r=radius;
  r2=r*r;
  diameter=r*2;
  volume=2*r2*acos(-1.)*len;
}
void cylinder::randompos(point &p) {
  double l=r2+1;
  p.z = slp.random_one()*len;
  while (l>r2) {
    p.x = slp.random_half()*diameter;
    p.y = slp.random_half()*diameter;
    l=p.x*p.x+p.y*p.y;
  };
}
string cylinder::info() {
  ostringstream o;
  o << container::info()
    << "#   Shape                = Cylindrical" << endl
    << "#   Lenght               = " << len <<endl
    << "#   Radius               = " << r << endl;
  return o.str();
}
string cylinder::povray() {
  ostringstream o;
  o << "cylinder {<0,0,0>,<0,0" << len << ">," << r <<" texture {cell}}\n"
    << "cylinder {<0,0,0>,<0,0" << len << ">,0.5 texture {cell}}\n";
  return o.str();
}     
}//namespace
