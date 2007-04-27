#include "potentials.h"

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
  int a,b,len,tmp;
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

/*********************************
   J'TH WITH ALL OTHER PARTICLES
 *********************************/
bool hardsphere::overlap(vector<particle> &p, int j) {
  int ps=p.size();
  for (int i=0; i<j; i++)
    if (p[i].overlap(p[j])==true) return true;
  for (int i=j+1; i<ps; i++)
    if (p[i].overlap(p[j])==true) return true;
  return false;
}
/********************************
   GROUP WITH REST OF PARTICLES
 ********************************/
bool hardsphere::overlap(vector<particle> &p, group &g) {
  short int n=g.beg, psize=p.size();
  for (short int i=0; i<n; i++)
    if (overlap(p, g, i)==true) return true;
  for (short int i=n+1; i<psize; i++) 
    if (overlap(p, g, i)==true) return true;
  return false;
}
//group with j'th particle
bool hardsphere::overlap(vector<particle> &p, group &g, int j) {
  if (g.beg==-1) return false;
  if (g.radius>0)
    if ( abs(g.cm.dist(p[j])) > g.radius+p[j].radius)
      return false;
  int len=g.end+1;

  if (g.find(j)==false) {          //check if j is part of g
    for (int i=g.beg; i<len; i++)
      if (p[i].overlap(p[j])==true)
        return true;
  } else {                              //it is! avoid self-overlap...
    for (int i=g.beg; i<j; i++)
      if (p[i].overlap(p[j])==true) return true;
    for (int i=j+1; i<len; i++)
      if (p[i].overlap(p[j])==true) return true;
  };
  return false;
}

/********************************
 *
 * GROUP <-> ARBITRARY PARTICLE
 *
 ********************************/
bool hardsphere::overlap(vector<particle> &p, group &g, particle &a) {
  int len=g.end+1;
  if (g.beg!=-1)
    for (int i=g.beg; i<len; i++)
      if (p[i].overlap(a)==true)
	return true;
  return false;
}


bool hardsphere::overlap(vector<particle> &p, group &g1, group &g2) {
  if (g1.beg==-1 || g2.beg==-1) return false;
  if (g1.radius!=0 && g2.radius!=0)
    if (abs(g1.cm.dist(g2.cm)) > g1.radius+g2.radius)
      return false;
  int ilen=g1.end+1, jlen=g2.end+1;
  for (int i=g1.beg; i<ilen; i++) {
    for (int j=g2.beg; j<jlen; j++) {
      if ( p[i].overlap(p[j])==true )
        return true;
    };
  };
  return false;
}


