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
  cout << "# WARNING! FILE " << file << " NOT READ!\n";
  return false;
}
bool io::writefile(string file, string s) {
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

//--------------- IOAAM ---------------------
ioaam::ioaam(species &spc) : iopart(spc) {}
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
vector<particle> ioaam::load(string file) {
  v.resize(0);
  p.resize(0);
  if (readfile(file,v)==true) {
    strip(v,"#");
    unsigned short n=atoi(v[0].c_str());
    for (unsigned short i=1; i<=n; i++)
      p.push_back( s2p(v[i]) );
  }
  return p;
}
bool ioaam::save(string file, vector<particle> &) {
  ostringstream o;
  o << p.size() << endl;
  for (unsigned short i=0; i<p.size(); i++)
    o << p2s(p[i]);
  return writefile(file, o.str());
}

//----------------- IOXYZ ----------------------
//ioxyz::ioxyz(species &spc) : iopart(spc) {
//}
ioxyz::ioxyz(species &s) : iopart(s) {
}
particle ioxyz::s2p(string &s) {
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

bool ioxyz::save(string file, vector<particle> &p) {
  ostringstream o;
  o << p.size() << endl << endl;
  for (unsigned short i=0; i<p.size(); i++)
    o << spcPtr->d[p[i].id].name << " "
      << p[i].x << " " << p[i].y << " " << p[i].z
      << endl;
  o << endl;
  return writefile(file, o.str());
}

vector<particle> ioxyz::load(string file) {
  v.resize(0);
  particle dummy;
  
  if (readfile(file,v)==true) {
    strip(v,"#");
    unsigned short n=atoi(v[0].c_str());
    if (n!=sys->p.size()) {
      cout << "!! System vector and .coord vector are out of sync! !!\n"
           << "!! Closing the factory since n=" <<n<< "and p.size() ="<<sys->p.size()<<" !!\n"
           << endl;
    }else{ 
    for (unsigned short i=1; i<=n; i++) {
      sys->p[i-1]=s2p(v[i]);
      sys->p[i-1].charge=s2p(v[i]).charge;
      sys->p[i-1].radius=s2p(v[i]).radius;
      sys->p[i-1].mw=s2p(v[i]).mw;
      }
    cout << "# Old coordinates loaded\n";
    }
  }else{ cout << "# Can't find .coord${jobid} \n";}
  return p;
}

//----------------- IOPOV ----------------------
//iopov::iopov(species &spc) : iopart(spc) {
iopov::iopov(container &c) : iopart(c) {
  clear();
  header();
  light(100.);
  camera();
  o << c.povray();
}

vector<particle> iopov::load(string s) {
  vector<particle> dummy(0);
  return dummy;
}

void iopov::clear() { o.str(""); }
void iopov::header() {
  clear();
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
    "}\n"
    "#declare cell=texture {\n"
    " pigment {color rgbf <1,1,1,.9>}\n"
    " finish {phong .1 ambient .2}\n"
    "}\n";
}
void iopov::camera() {
  o << "camera {" << endl
    << "  location <100,100,0>" << endl
    << "  look_at <0,0,0>\n}" << endl;
}
void iopov::light(float r) {
  o << "light_source {\n"
    << "  <"<<r*1.5<<","<<r*1.5<<","<<r*0.5<<">\n"
    << "  color rgb <1,1,1>\n}\n";
  o << "light_source {\n"
    << "  <"<<-r*1.5<<","<<r*1.5<<","<<-r*0.5<<">\n"
    << "  color rgb <1,1,1>\n}" << endl;
}
string iopov::p2s(particle &p, int i) {
  string tex;
  ostringstream s;
  if (p.charge>0)  tex="redish";
  if (p.charge<0)  tex="greyish";
  if (p.charge==0) tex="white";
  if (p.radius>0 && p.id!=particle::GHOST) 
    s << "sphere {<"<<p.x<<","<<p.y<<","<<p.z<<">,"<<p.radius
      << " texture {"<<tex<<"}}\n";
  return s.str();
}
bool iopov::save(string file, vector<particle> &p) {
  for (unsigned short i=0; i<p.size(); i++)
    o << p2s(p[i]);
  return writefile(file, o.str());
}
void iopov::connect(point &p1, point &p2, float radius) {
  o << "cylinder {\n"
    << "  <"<<p1.x<<","<<p1.y<<","<<p1.z<<">,"
    << "  <"<<p2.x<<","<<p2.y<<","<<p2.z<<">,"<<radius
    << "  texture {redish}\n }\n";
}


/*PQR
// 1234567890123456789012345678901234567890123456789012345678901234567890
// ATOM     25 H1A  CR      4      23.215  17.973  15.540 -0.017 1.000
*/

//----------------- IOXTC ----------------------
#ifdef GROMACS
ioxtc::ioxtc(container::container &con, float len) : iopart(con) {
  time=step=0;
  setbox(len);
  xd=open_xtc("coord.xtc", "w");
}

void ioxtc::setbox(float len) {
  for (char i=0; i<3; i++)
    for (char j=0; j<3; j++)
      box[i][j]=0;
  box[0][0]=len/10.; // corners of the
  box[1][1]=len/10.; // rectangular box
  box[2][2]=len/10.; // (in nanometers!)
};

bool ioxtc::save(string file, vector<particle> &p) {
  if (p.size()<3300) {
    for (unsigned short i=0; i<p.size(); i++) {
      x[i][0]=p[i].x/10;      // AA->nm
      x[i][1]=p[i].y/10;
      x[i][2]=p[i].z/10;
    }
    write_xtc(xd,p.size(),step++,time++,box,x,3300.);
  }
}
void ioxtc::close() { close_xtc(xd); }
#endif


