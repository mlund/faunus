#include "potentials.h"

pot_setup::pot_setup() {
  lB=7.1;
  eps=2;
  epsi=2;
  epso=80;
}

pot_setup::pot_setup(inputfile &in) {
  lB    = in.getflt("bjerrum", 7.1);
  eps   = in.getflt("LJeps", 2);
  box   = in.getflt("boxlen");
  kappa = in.getflt("kappa");
  epsi  = in.getflt("epsi",2);
  epso  = in.getflt("epso",80);
  a     = in.getflt("cavity");
  hydroscale = in.getflt("hydroscale", 4.0);
}

string pot_minimage::info() {
  ostringstream o;
  o << pot_lj::info()
    << "#   Bjerrum length    = " << f << endl
    << "#   Image length      = " << box << endl;
  return o.str();
}

string pot_hscoulomb::info() {
  ostringstream o;
  o << pot_hs::info()
    << "#   Bjerrum length    = " << f << endl;
  return o.str();
}
string pot_coulomb::info() {
  ostringstream o;
  o << pot_lj::info()
    << "#   Bjerrum length    = " << f << endl;
  return o.str();
}

string pot_debyehuckel::info() {
  ostringstream o;
  o << pot_lj::info()
    << "#   Bjerrum length    = " << f     << endl
    << "#   Debye length      = " << 1./k  << endl;
  return o.str();
}

string pot_debyehuckelP3::info() {
  ostringstream o;
  o << pot_lj::info()
    << "#   Bjerrum length    = " << f     << endl
    << "#   Debye length      = " << 1./k  << endl;
  return o.str();
}

/*!\param spc Species class.
 * \param pmfir Directory in which to search for PMF's
 * \param type Particle name to search for. I.e. "NA" or "CL"
 */
void pot_datapmf::loadpmf(species &spc, string pmfdir, string type) {
  string n1,n2;                                                      
  unsigned short i,j=spc.id(type); 
  for (i=particle::FIRST; i<particle::LAST; i++) { 
    n1=spc.d[i].name+"-"+spc.d[j].name;                                      
    n2=spc.d[j].name+"-"+spc.d[i].name;                                      
    if (loadpmf( spc, pmfdir + n1 + ".dat")==false)                       
      loadpmf( spc, pmfdir + n2 + ".dat");                                
  };                                   
}

// Show table of loaded PMF's
void pot_datapmf::showpmf(species &spc) {
  cout << "# --- LOADED PMF's ----------------------------------\n"; 
  cout << "# (a,b,resolution,r_max,file)\n"; 
  for (unsigned int i=particle::FIRST; i<particle::LAST; i++)
    for (unsigned int j=particle::FIRST; j<particle::LAST; j++)
      if (pmfd[i][j].xmax>0) 
        cout << "# " << spc.d[i].name << " " << spc.d[j].name << " "
          << pmfd[i][j].res << " "
          << pmfd[i][j].xmax << " "
          << pmfd[i][j].comment << endl; 
  cout << endl; 
}

/*! Load PMF(s) from a file. File format: Each set starts
 * with "#$ type1 type2 length". Several sets can be present
 * in the same file.
 */
bool pot_datapmf::loadpmf(species &spc, string filename) {
  string s,a_str,b_str;
  int a,b,len;
  vector<double> x,y;
  ifstream fh(filename.c_str());
  if (fh) {
    //scan file, word-by-word
    while (!fh.eof()) {
      fh >> s;
      if (s.find("#$")==0) {
        s.clear();
        fh >> a_str >> b_str >> len;
        x.resize(len);
        y.resize(len);
        a = spc.id(a_str);
        b = spc.id(b_str);
        if (a>b)
          swap(a,b);
        for (int i=0; i<len; i++) {
          fh >> x[i] >> y[i];
          y[i]=y[i]/f; // unit: kT/f
        };
        if (pmfd[a][b].xmax==0) {
          pmfd[a][b].add( x, y );
          pmfd[a][b].comment=filename;
        };
      };
    };
    fh.close();
    return true;
  }
  return false;
}


