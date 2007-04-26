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
double interaction::energy(vector<particle> &p, int j) {
  int ps=p.size();
  double u=0;
  for (int i=0; i<j; i++)
    u+=pairpot( p[i],p[j] );
  for (int i=j+1; i<ps; i++)
    u+=pairpot( p[i],p[j] );
  return f*u;
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
double interaction::energy(vector<particle> &p, group &g) {
  int n=g.end+1, psize=p.size();
  double u=0;
  for (int i=g.beg; i<n; i++) {
    for (int j=0; j<g.beg; j++)
      u += pairpot(p[i],p[j]);
    for (int j=n; j<psize; j++)
      u += pairpot(p[i],p[j]);
  };
  return f*u;
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

double interaction::energy(vector<particle> &p, group &g, int j) {
  double u=0;
  int len=g.end+1;
  if (g.find(j)==true) {   //avoid self-interaction...
    for (int i=g.beg; i<j; i++)
      u+=pairpot(p[i],p[j]);
    for (int i=j+1; i<len; i++)
      u+=pairpot(p[i],p[j]);
  } else                        //simple - j not in g
    for (int i=g.beg; i<len; i++)
      u+=pairpot(p[i],p[j]);
  return f*u;  
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

double interaction::energy(vector<particle> &p, group &g, particle &a) {
  if (g.beg==-1)
    return 0;
  double u=0;
  int len=g.end+1;
  for (int i=g.beg; i<len; i++) 
    u+=pairpot(a, p[i]); 
  return f*u;
}

/*********************
   SYSTEM ENERGY
 *********************/
double interaction::energy(vector<particle> &p) {
  double u=0;
  int n = p.size();
  for (int i=0; i<n-1; i++)
    for (int j=i+1; j<n; j++)
      u += pairpot(p[i], p[j]);
  return f*u; 
}


/**********************
   BETWEEN TWO GROUPS
 **********************/
double interaction::energy(vector<particle> &p, group &g1, group &g2) {
  int ilen=g1.end+1; 
  int jlen=g2.end+1;
  double u=0;
  for (int i=g1.beg; i<ilen; i++)
    for (int j=g2.beg; j<jlen; j++)
      u += pairpot(p[i],p[j]);
  return f*u;
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

//non-electrostatic energy of chain particle i with rest of the chain
double interaction::chain(vector<particle> &p, group &g, int i) {
  double u=0;
  //the first ?
  if (i==g.beg) {
    u+=quadratic( p[i], p[i+1] );
    if (g.graftpoint>-1)
      u+=quadratic( p[i], p[g.graftpoint] );
    return u;
  };

  //the last ?
  if (i==g.end)
    return quadratic( p[i], p[i-1] );

  //otherwise...
  return quadratic(p[i], p[i+1]) + quadratic(p[i], p[i-1]);
}

//graft energy of the chain, g.
double interaction::graft(vector<particle> &p, group &g) {
  if (g.graftpoint!=-1)
    return quadratic( p[g.graftpoint], p[g.beg]);
  else
    return 0;
}

/*!
 * ...between the two dipoles a and b, separated by the
 * distance r.
 * \f$ \beta u(r) = l_B \frac{a_x b_x + a_y b_y - 2a_z b_z  }{r^3}\f$
 */
double interaction::dipdip(point &a, point &b, double r) {
  return f*( a.x*b.x + a.y*b.y - 2*a.z*b.z )/(r*r*r);
}
double interaction::iondip(point &a, double q, double r) {
  return -f*q*a.z/(r*r);
}

// Total electrostatic potential in a point
double interaction::pot(vector<particle> &p, point &a) {
  double u=0;
  int ps=p.size();  
  for (int i=0; i<ps; i++)
    u+=p[i].charge / p[i].dist(a);
  return f*u;
}

// Internal (NON)-electrostatic energy in group
double interaction::internal(vector<particle> &p, group &g) {
  int glen=g.end+1;
  double u=0;
  if (g.beg==-1)
    return 0;
  
  if (g.chain==true) {
    for (int i=g.beg; i<g.end; i++)
      u+=quadratic( p[i], p[i+1]);
    if (g.graftpoint>-1)
      u+=quadratic( p[g.beg], p[g.graftpoint] );
    return u;
  }
  else {
    for (int i=g.beg; i<glen-1; i++)
      for (int j=i+1; j<glen; j++)
        u+=pairpot(p[i],p[j]);
  }
  return f*u;
}

