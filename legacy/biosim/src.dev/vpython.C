#include "vpython.h"

vpython::vpython() {
  ostr << "#!/sw/bin/python\n\n"
       << "from visual import *\n\n";
};

void vpython::add(vector<particle> &p, group &g) {
  string col;
  for (int i=g.beg; i<=g.end; i++)
    if (p[i].radius>0) {
      if (p[i].charge==0) col="color.white";
      if (p[i].charge<0)  col="color.blue";
      if (p[i].charge>0) col="color.red";
      ostr << setiosflags( ios::fixed );
      ostr.precision(1);
      ostr << "sphere(pos=(" << p[i].x << "," << p[i].y << "," << p[i].z << "),"
           << "radius=" << p[i].radius << ","
           << "color=" << col << ")\n";
    };
};

bool vpython::save(string filename) {
  bool rc=false;
  ofstream f( filename.c_str() );
  if (f) {
    f << ostr.str();
    f.close();
    rc=true;
  };
  return rc;
};


////////

vrml::vrml() {
  ostr << "#!/sw/bin/python\n"
       << "from Scientific.Visualization.VRML2 import *\n\n"
       << "s=Scene([])\n\n"
       << "scale = ColorScale(10.)\n"
       << "neuCol=scale(1.)\n"
       << "neu=Material(diffuse_color=neuCol,shininess=0.9)\n"
       << "catCol=scale(2.)\n"
       << "cat=Material(diffuse_color=catCol,shininess=0.9)\n"
       << "anCol=scale(3.)\n"
       << "an=Material(diffuse_color=anCol,shininess=0.9)\n";
};

void vrml::add(vector<particle> &p, group &g) {
  string col;              
  for (int i=g.beg; i<=g.end; i++)
    if (p[i].radius>0) {   
      if (p[i].charge==0) col="neu";
      if (p[i].charge<0)  col="an";
      if (p[i].charge>0) col="cat";
      ostr << setiosflags( ios::fixed );
      ostr.precision(1);
      ostr << "s.addObject(Sphere(Vector("
           << p[i].x << "," << p[i].y << "," << p[i].z
           << ")," << p[i].radius << ",material="<< col<< "))\n";
    };
};

bool vrml::save(string filename) {
  bool rc=false;
  ofstream f( filename.c_str() );
  if (f) {
    ostr << "s.writeToFile(\"test.wrl\")" << endl;
    f << ostr.str();
    f.close();
    rc=true;
  };
  return rc;
};

