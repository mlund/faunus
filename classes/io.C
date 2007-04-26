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

vector<particle> iofile::load(string file) {
}

ioaam::ioaam(species &spc) {
  spcPtr = &spc;
}

particle ioaam::s2p(string &s) {
  stringstream o;
  particle p;
  string name, num;
  o << s;
  o >> name >> num
    >> p.x >> p.y >> p.z
    >> p.charge >> p.mw >> p.radius;
  p.id = spcPtr->id(name); 
  return p;
}

/*
// 1234567890123456789012345678901234567890123456789012345678901234567890
// ATOM     25 H1A  CR      4      23.215  17.973  15.540 -0.017 1.000

void io::pqrline(particle &) {
};

void io::savepqr(string filename, vector<particle> &p) {
};
*/
