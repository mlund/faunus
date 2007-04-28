#include "io.h"

/*!
 * The first line in the AAM file format is the number of
 * particles. The following lines defines the particles according
 * to: \n
 *   name num x y z charge weight radius
 */

//-------------- IO -------------------
bool io::readfile(string file, vector<string> &v) {
  string s;
  ifstream f(file.c_str() );
  if (f) {
    while (getline(f,s)) 
      v.push_back(s);
    f.close();
    return true;
  }
  return false;
}

//--------------IO PARTICLE ------------------
iopart::iopart(species &spc) { spcPtr=&spc; }

vector<particle> iopart::load(string file) {
  unsigned int beg,end;
  vector<string> v;
  vector<particle> f,p;
  if (readfile(file,v)==true) {
    beg=first();
    end=last();
    for (int i=beg; i<=end; i++)
      p.push_back( s2p(v[i]) );
    return p;
  }
}

bool iopart::save(vector<particle> &p, string file) {
}

//--------------- IOAAM ---------------------
ioaam::ioaam(species &spc) : iopart(spc) {}
unsigned int ioaam::first() { return 1; }
unsigned int ioaam::last() { return atoi(v[0].c_str() ); }

particle ioaam::s2p(string &s) {
  particle p;
  stringstream o;
  string name, num;
  o << s;
  o >> name >> num
    >> p.x >> p.y >> p.z
    >> p.charge >> p.mw >> p.radius;
  p.id = spcPtr->id(name); 
  return p;
}

string ioaam::p2s(particle &) {}

//----------------- IOPOV ----------------------
string iopov::p2s(particle &p) {
  stringstream s;
  string tex;
  if (p.charge>0)  tex="redish";
  if (p.charge<0)  tex="greyish";
  if (p.charge==0) tex="white";
  if (p.radius>0 && p.id!=particle::GHOST) 
    s << " sphere {<"<<p.x<<","<<p.y<<","<<p.z<<">,"<<p.radius
      << " texture {"<<tex<<"}}\n";
  return s.str();
}

string iopov::header() {
  string s(
    "#declare white=texture {\n"
    " pigment {color rgb <1,1,1>}\n"
    " finish {phong .9 ambient .1 reflection 0.2}\n"
    "}\n"
    "#declare black=texture {\n"
    " pigment {color rgb <0,0,0>}\n"
    " finish {phong .9 ambient .2 reflection .2}\n"
    "}\n"
    "#declare transp=texture {\n"
    " pigment {color rgbf <1,1,1,.9>}\n"
    " finish {phong .1 ambient .2}\n"
    "}\n"
    "#declare greyish=texture {\n"
    " pigment {color rgbf <.5,.5,.5,.7>}\n"
    " finish {phong .9 ambient .1 reflection .2}\n"
    "}\n"
    "#declare redish=texture {\n"
    " pigment {color rgb <1,0,0>}\n"
    " finish {phong .9 ambient .1 reflection .2}\n"
    "}\n" );
  return s;
}

/*
// 1234567890123456789012345678901234567890123456789012345678901234567890
// ATOM     25 H1A  CR      4      23.215  17.973  15.540 -0.017 1.000

void io::pqrline(particle &) {
};

void io::savepqr(string filename, vector<particle> &p) {
};
*/
