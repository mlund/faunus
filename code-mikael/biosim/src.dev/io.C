#include "io.h"

// net charge
double io::charge(vector<particle> &p) {
  double q=0;
  for (unsigned int i=0; i<p.size(); i++)
    q+=p[i].charge;
  return q;
};

// center of mass
point io::cm(vector<particle> &p) {
  particle cm;
  double x=0,y=0,z=0,m=0;
  for (int i=0; i<p.size(); i++) {
    x+=p[i].x*p[i].mw;
    y+=p[i].y*p[i].mw;
    z+=p[i].z*p[i].mw;
    m+=p[i].mw;
  };
  cm.x = x/m;
  cm.y = y/m;
  cm.z = z/m;
  return cm;
};

// load file w. space separated fields formatted as
// name num x y z charge weight radius
vector<particle> io::loadaam(species &spc, string filename) {
  vector<particle> p(0);
  point c;
  int n, num;
  string s,res;
  ifstream f(filename.c_str());
  if (f) {
    do getline(f,s);
    while (s.find("#")!=string::npos);
    n=atoi(s.c_str());
    p.resize(n);
    c.y=-6.*n/2.;
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
};

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
};

/*
// 1234567890123456789012345678901234567890123456789012345678901234567890
// ATOM     25 H1A  CR      4      23.215  17.973  15.540 -0.017 1.000

void io::pqrline(particle &) {
};

void io::savepqr(string filename, vector<particle> &p) {
};
*/
