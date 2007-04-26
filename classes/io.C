#include "io.h"

/*!
 * The first line in the AAM file format is the number of
 * particles. The following lines defines the particles according
 * to: \n
 *   name num x y z charge weight radius
 */
vector<particle> io::loadaam(species &spc, string filename) {
  vector<particle> p(0);
  int n, num;
  string s,res;
  ifstream f(filename.c_str());
  if (f) {
    do getline(f,s);
    while (s.find("#")!=string::npos);
    n=atoi(s.c_str());
    p.resize(n);
    for (int i=0; i<n; i++) {
      f >> res >> num >> p[i].x >> p[i].y
        >> p[i].z >> p[i].charge >> p[i].mw >> p[i].radius;
      p[i].id = spc.id(res);
    };
    f.close();
  }
  else
    cout << "# Error! File '"<<filename<<"' not found.";
  return p;
}

// save group
bool io::saveaam(species &spc, string file,
                      vector<particle> &p, group &g) {
  ofstream f( file.c_str() );
  if (f) {
    f << g.size() << endl;
    for (int i=g.beg; i<=g.end; i++) {
      if (p[i].radius>0) {
        if (p[i].id<spc.d.size())
          f << spc.d[p[i].id].name <<" ";
        else 
          f << "UNK" << " ";
        f << i+1 << " ";
        f << p[i].x << " " << p[i].y << " " << p[i].z << " ";
        f << p[i].charge << " " << p[i].mw << " " << p[i].radius << endl;
      };
    };
    f.close();
    return true;
  };
  return false;
}

//-------------- IO FILE -------------------
iofile::iofile(species &spc) { spcPtr=&spc; }

bool iofile::readfile(string file) {
  string s;
  v.resize(0);
  ifstream f(file.c_str() );
  if (f) {
    while (!f.eof()) {
      getline(f,s);
      if (s.size()>0 && s.find("#")==string::npos)
        v.push_back(s);
    }
    f.close();
    return true;
  }
  return false;
}

vector<particle> iofile::load(string file) {
  unsigned int beg,end;
  vector<particle> p(0);
  if (readfile(file)==true) {
    beg=first();
    end=last();
    for (int i=beg; i<=end; i++)
      p.push_back( s2p(v[i]) );
    return p;
  }
}

//--------------- IOAAM ---------------------
ioaam::ioaam(species &spc) : iofile(spc) {}
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
