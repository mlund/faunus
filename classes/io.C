#include "io.h"

/*!
 * The first line in the AAM file format is the number of
 * particles. The following lines defines the particles according
 * to: \n
 *   name num x y z charge weight radius
 */

//------------------ IO -------------------
bool io::readfile(string file, vector<string> &v) {
  ifstream f(file.c_str() );
  if (f) {
    while (getline(f,s)) 
      v.push_back(s);
    f.close();
    return true;
  }
  return false;
}
bool io::writefile(string file, string &s) {
  ofstream f(file.c_str());
  if (f) {
    f << s;
    f.close();
    return true;
  }
  return false;
}
void io::strip(vector<string> &v, string pat) {
  vector<string>::iterator iter=v.begin();
  while (iter!=v.end())
    if ((*iter).find(pat)!=string::npos)
      v.erase(iter);
    else iter++;
}

//--------------IO PARTICLE ------------------
iopart::iopart(species &spc) { spcPtr=&spc; }
vector<particle> iopart::load(string file) {
  v.resize(0);
  p.resize(0);
  if (readfile(file,v)==true) {
    strip(v, "#");
    for (unsigned char i=first(); i<=last(); i++)
      p.push_back( s2p(v[i]) );
  }
  return p;
}
bool iopart::save(vector<particle> &p, string file) {
  string s=header(p);
  for (unsigned char i=0; i<p.size(); i++)
    s+=p2s(p[i], i);
  return writefile(file, s);
}

//--------------- IOAAM ---------------------
ioaam::ioaam(species &spc) : iopart(spc) {}
unsigned char ioaam::first() { return 1; }
unsigned char ioaam::last() { return atoi(v[0].c_str())+first()-1; }
string ioaam::p2s(particle &p, int i) {
  ostringstream o;
  o.precision(30);
  o << spcPtr->d[p.id].name<<" "<<i+1<<" "
    << p.x<<" "<< p.y<<" "<<p.z<<" "
    << p.charge<<" "<<p.radius<<endl;
  return o.str();
}
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
string ioaam::header(vector<particle> &p) {
  ostringstream o;
  o << p.size() << endl;
  return o.str();
}

//----------------- IOPOV ----------------------
iopov::iopov(species &spc) : iopart(spc) {
  o << 
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
    "}\n";
}
string iopov::header(vector<particle> &p, int i) { return o.str(); }
string iopov::p2s(particle &p, int i) {
  string tex;
  ostringstream s;
  if (p.charge>0)  tex="redish";
  if (p.charge<0)  tex="greyish";
  if (p.charge==0) tex="white";
  if (p.radius>0 && p.id!=particle::GHOST) 
    s << " sphere {<"<<p.x<<","<<p.y<<","<<p.z<<">,"<<p.radius
      << " texture {"<<tex<<"}}\n";
  return s.str();
}


/*
// 1234567890123456789012345678901234567890123456789012345678901234567890
// ATOM     25 H1A  CR      4      23.215  17.973  15.540 -0.017 1.000

void io::pqrline(particle &) {
};

void io::savepqr(string filename, vector<particle> &p) {
};
*/
