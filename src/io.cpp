#include "faunus/io.h"
namespace Faunus {
  /*!
   * The first line in the AAM file format is the number of
   * particles. The following lines defines the particles according
   * to: \n
   *   name num x y z charge weight radius
   */

  //------------------ IO -------------------
  bool io::readfile(string file, vector<string> &v) {
    std::ifstream f(file.c_str() );
    if (f) {
      while (getline(f,s))
        v.push_back(s);
      f.close();
      return true;
    }
    std::cerr << "# WARNING! FILE " << file << " NOT READ!\n";
    return false;
  }
  bool io::writefile(string file, string s) {
    std::ofstream f(file.c_str());
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
  void io::splash(string f) {
    vector<string> t;
    t.clear();
    readfile(f,t);
    for (int i=0; i<t.size(); i++)
      std::cout << "# "<<t[i]<<endl;
  }

  /*
     bool io_aam::load(container &c, string file) {
     return false;
     }
     */

  //--------------- IOAAM ---------------------
  ioaam::ioaam() {}
  string ioaam::p2s(particle &p, int i) {
    std::ostringstream o;
    o.precision(5);
    o << atom[p.id].name << " " << i+1 << " "
      << p.x<<" "<< p.y<<" "<<p.z<<" "
      << p.charge << " " << p.mw << " " << p.radius << std::endl;
    return o.str();
  }
  particle ioaam::s2p(string &s) {
    particle p;
    std::stringstream o;
    string name, num;
    o << s;
    o >> name >> num
      >> p.x >> p.y >> p.z
      >> p.charge >> p.mw >> p.radius;
    p.id = atom[name].id;
    p.hydrophobic = atom[p.id].hydrophobic;
    return p;
  }
  vector<particle> ioaam::load(string file) {
    vector<string> v;
    vector<particle> p;
    if (readfile(file,v)==true) {
      strip(v,"#");
      unsigned int n=atoi(v[0].c_str());
      for (unsigned int i=1; i<=n; i++)
        p.push_back( s2p(v[i]) );
    }
    v.clear();
    return p;
  }

  bool ioaam::load(container &con, string file) {
    vector<particle> a=load(file);
    if (a.size()==con.p.size()) {
      con.p = a;
      con.trial = a;
      std::cout << "# " << a.size() << " particles loaded from " << file << endl;
      return true;
    }
    std::cout << "# " << file << " skipped due to container mismatch!\n";
    return false;
  }

  bool ioaam::save(string file, vector<particle> &p) {
    std::ostringstream o;
    o << p.size() << std::endl;
    for (unsigned int i=0; i<p.size(); i++)
      o << p2s(p[i], i);
    return writefile(file, o.str());
  }

  /*!
   * Searches inputfile object for the following keywords:\n
   *   "nprot#" -- number protein # structures\n
   *   "protein#" -- name of protein # structure file\n
   * The found proteins will be inserted in the particle vector
   * and for each structure a group will be appended to the
   * macromolecular vector.
   *
   * \author Mikael Lund
   */
  void ioaam::load(container &con, inputfile &in, vector<macromolecule> &g) {
    short cnt=1,nprot;
    do {
      std::ostringstream os_prot,os_file;
      os_prot << "nprot" << cnt;
      os_file << "protein" << cnt++;
      nprot = in.getint(os_prot.str(), 0);
      if (nprot>0)
        for (short i=0; i<nprot; i++) {
          macromolecule m;
          m.add(con, load( in.getstr(os_file.str() )  ), true ); 
          m.name = in.getstr(os_file.str());
          m.masscenter(con);
          g.push_back(m);
        }
    } while (nprot>0);
  }
  void ioaam::loadlattice(container &con, inputfile &in, vector<macromolecule> &g) {
    std::ostringstream n_prot;
    short pcnt=1, n, N=0, unitclen;
    double len;
    n_prot << "nprot" <<pcnt++;
    n=in.getint(n_prot.str(),0);  //Determine the number of proteins
    if (n!=0) {
      N+=n;
      n_prot.flush();
      n_prot << "nprot" <<pcnt++;
      n=in.getint(n_prot.str(),0);
    }
    len=pow(con.getvolume(),1./3);    //Calculate the length of a unit cell
    unitclen=int( pow(N,1./3) );
    len/=(unitclen+1);
    short cnt=1,nprot;
    do {
      std::ostringstream os_prot,os_file;
      os_prot << "nprot" << cnt;
      os_file << "protein" << cnt++;
      nprot = in.getint(os_prot.str(), 0);
      if (nprot>0)
        for (short i=0; i<nprot; i++) {
          macromolecule m;
          m.add(con, load(in.getstr(os_file.str())), false ); 
          m.masscenter(con);
          m.center(con);
          m.name = in.getstr( os_file.str() );
          g.push_back(m);
        }
    } while (nprot>0);
    point P;
    n=0;
    for (short p=0;p<unitclen+1;p++){
      P.x=len*p;
      for (short q=0; q<unitclen+1;q++){
        P.y=len*q;
        for (short r=0; r<unitclen+1;r++){
          P.z=len*r;
          if (n<g.size()) {
            g[n].move(con, P);
            g[n].accept(con);
            n++;
          }
        }
      }
    }
  }

  //----------------- IOXYZ ----------------------
  ioxyz::ioxyz() {}
  particle ioxyz::s2p(string &s) {
    particle p;
    std::stringstream o;
    string name, num;
    o << s;
    o >> name >> num
      >> p.x >> p.y >> p.z
      >> p.charge >> p.mw >> p.radius;
    p.id = atom[name].id;
    p.hydrophobic = atom[p.id].hydrophobic;
    return p;
  }

  bool ioxyz::save(string file, vector<particle> &p) {
    std::ostringstream o;
    o << p.size() << std::endl << std::endl;
    for (unsigned short i=0; i<p.size(); i++)
      o << atom[p[i].id].name << " "
        << p[i].x << " " << p[i].y << " " << p[i].z
        << std::endl;
    o << std::endl;
    return writefile(file, o.str());
  }

  vector<particle> ioxyz::load(string file) {
    v.resize(0);
    particle dummy;

    if (readfile(file,v)==true) {
      strip(v,"#");
      unsigned short n=atoi(v[0].c_str());
      if (n!=sys->p.size()) {
        std::cerr << "!! System vector and .coord vector are out of sync! !!\n"
          << "!! Closing the factory since n=" <<n<< "and p.size() ="<<sys->p.size()<<" !!\n"
          << std::endl;
      }else{ 
        for (unsigned short i=1; i<=n; i++) {
          sys->p[i-1]=s2p(v[i]);
          sys->p[i-1].charge=s2p(v[i]).charge;
          sys->p[i-1].radius=s2p(v[i]).radius;
          sys->p[i-1].mw=s2p(v[i]).mw;
        }
        std::cout << "# Old coordinates loaded\n";
      }
    }else{ std::cout << "# Can't find .coord${jobid} \n";}
    return p;
  }
  //----------------- IOPQR ----------------------
  /*!
   *
   *
   */
  iopqr::iopqr() { }
  bool iopqr::save(string file, vector<particle> &p) {
    string name;
    int i, nres=1, natom=1;
    char buf[100];
    std::ostringstream o;
    for (i=0; i<p.size(); i++) {
      // index, atom->name, atom->resname, atom->resid,x, y, z, atom->charge, atom->radius
      name=atom[p[i].id].name;
      sprintf(buf, "ATOM  %5d %-4s %s %5d    %8.3f %8.3f %8.3f %.3f %.3f\n",
          natom++, name.c_str(), name.c_str(), nres,
          p[i].x, p[i].y, p[i].z, p[i].charge, p[i].radius );
      o << buf;
      if ( atom[p[i].id].name=="CTR" ) nres++;
    }
    return writefile(file, o.str());
  }
  bool iopqr::save(string file, vector<particle> &p, titrate &tit) {
    vector<particle> tmp = p;
    for (unsigned short i=0; i<p.size(); i++)
      tmp[i].charge = tit.avgcharge(p, i);
    return save(file, tmp);
  }

  //----------------- IOGRO ----------------------
  /*!
   * \param s Species class for particle recognition
   * \param in Inputfile. Keyword "boxlen" is used to specify a cubic box
   */
  iogro::iogro(inputfile &in) {
    len = in.getflt("boxlen", 0)/10;
  }
  bool iogro::save(string file, box &b) {
    len = b.len/10.;
    return save(file, b.p);
  }
  bool iogro::save(string file, vector<particle> &p) {
    string name;
    int i, nres=1, natom=1;
    char buf[79];
    std::ostringstream o;
    o << "# Generated by Faunus -- http://faunus.sourceforge.net"
      << std::endl << p.size() << std::endl;
    for (i=0; i<p.size(); i++) {
      name=atom[p[i].id].name;
      sprintf(buf, "%5d%5s%5s%5d%8.3f%8.3f%8.3f\n",
          nres,name.c_str(),name.c_str(),natom++,
          p[i].x/10,p[i].y/10,p[i].z/10);
      o << buf;
      if ( atom[p[i].id].name=="CTR" ) nres++;
    }
    if (len>0)
      o << len << " " << len << " " << len << std::endl;
    return writefile(file, o.str());
  }
  particle iogro::s2p(string &s) {
    std::stringstream o;
    string name;
    double x,y,z;
    o << s.substr(10,5) << s.substr(20,8) << s.substr(28,8) << s.substr(36,8);
    o >> name >> x >> y >> z;
    particle p = atom(name); 
    p.x=x*10; // nm->angstrom
    p.y=y*10;
    p.z=z*10;
    return p;
  }
  vector<particle> iogro::load(string file) {
    p.clear();
    v.resize(0);
    if (readfile(file,v)==true) {
      int last=atoi(v[1].c_str())+1;
      for (int i=2; i<=last; i++)
        p.push_back( s2p(v[i]) );
    }
    return p;
  }

  //----------------- IOPOV ----------------------
  iopov::iopov(container &c) {
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
    o << "camera {" << std::endl
      << "  location <100,100,0>" << std::endl
      << "  look_at <0,0,0>\n}" << std::endl;
  }
  void iopov::light(float r) {
    o << "light_source {\n"
      << "  <"<<r*1.5<<","<<r*1.5<<","<<r*0.5<<">\n"
      << "  color rgb <1,1,1>\n}\n";
    o << "light_source {\n"
      << "  <"<<-r*1.5<<","<<r*1.5<<","<<-r*0.5<<">\n"
      << "  color rgb <1,1,1>\n}" << std::endl;
  }
  string iopov::p2s(particle &p, int i) {
    string tex;
    std::ostringstream s;
    if (p.charge>0)  tex="redish";
    if (p.charge<0)  tex="greyish";
    if (p.charge==0) tex="white";
    if (p.radius>0) 
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
  ioxtc::ioxtc(float len) {
    time=step=0;
    setbox(len);
    xd=xdrfile_open("coord.xtc", "w");
  }
  vector<particle> ioxtc::load(string s) {
    vector<particle> dummy;
    return dummy;
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
    if (p.size()<5000) {
      for (unsigned short i=0; i<p.size(); i++) {
        x[i][0]=p[i].x/10;      // AA->nm
        x[i][1]=p[i].y/10;
        x[i][2]=p[i].z/10;
      }
      write_xtc(xd,p.size(),step++,time++,box,x,1300.);
    }
    return true;
  }
  void ioxtc::close() { xdrfile_close(xd); }

  bool ioxtc::OpenTrajectory(string s) {
    xd=xdrfile_open(s.c_str(), "r");
    return (xd==NULL) ? false : true;
  }
  /*
     bool ioxtc::LoadFrame(int i, vector<particle> &p, double boxlen &box) {
     int rc,natoms_xtc,step_xtc;
     float time_xtc;
     matrix box_xtc;
     rvec *x_xtc;
     float prec_xtc = 1000.0;
     while(1)
     {
  //rc = read_xtc(xd_read, natoms_xtc, &step_xtc, &time_xtc,box_xtc, x_xtc, &prec_xtc);
  //if (rc == 0) {
  //}
  //else break;
  }
  }
  */
};//namespace
