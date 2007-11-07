/*
 * Interact class.
 * (intermolecular interactions)
 *
 * M. Lund, 2004
 *
 */

#include "interact.h"

bool collision::celloverlap(vector<particle> &p, group &g, double cell_r) {
  double cell_r2 = cell_r*cell_r;  //radius squared
  int end=g.end+1;
  for (int i=g.beg; i<end; i++)
    if (p[i].x*p[i].x + p[i].y*p[i].y + p[i].z*p[i].z > cell_r2)
      return true;
  return false;
};

//check for overlapping charges within group - the threshold distance is
//set by "contact" (Aangstrom).
bool collision::chgoverlap(vector<particle> &p, group &g, double contact) {
  double s2=contact*contact;
  for (int i=g.beg; i<g.end; i++)
    for (int j=i+1; j<g.end+1; j++)
      if (p[i].charge!=0 && p[j].charge!=0)
        if (p[i].sqdist(p[j]) < s2)
	  return true;
  return false;
};

interact::interact(double bjerrum) {
  lB = bjerrum; // bjerrum length in aangstroms
};

// arbitrary particle with particle vector
bool collision::overlap(vector<particle> &p, particle &a) {
  int ps=p.size();
  for (int i=0; i<ps; i++)
    if (p[i].overlap(a)==true)
      return true;
  return false;
}; 

// internal collision of a subset of "p" defined by "g". Contact value (Aangstrom) taken
// from variable "contact".
bool collision::overlap(vector<particle> &p, vector<short int> &g, double contact) {
  double r2=contact*contact;
  int size=g.size();
  for (int i=0; i<size-1; i++)
    for (int j=i+1; j<size; j++)
      if (p[g[i]].sqdist(p[g[j]])<r2)
        return true;
  return false;
};

/**********************************
   PARTICLE <-> LIST OF PARTICLES
 **********************************/
double interact::energy(vector<particle> &p, int j, vector<short int> &v) {
  double u=0;
  int vs=v.size();
  for (int i=0; i<vs; i++)
    u+=energy(p[v[i]],p[j]);
  return lB*u;
};

/*********************************
   J'TH WITH ALL OTHER PARTICLES
 *********************************/
bool collision::overlap(vector<particle> &p, int j) {
  int ps=p.size();
  for (int i=0; i<j; i++)
    if (p[i].overlap(p[j])==true) return true;
  for (int i=j+1; i<ps; i++)
    if (p[i].overlap(p[j])==true) return true;
  return false;
};
double interact::energy(vector<particle> &p, int j) {
  if (p[j].charge==0) return 0;
  int ps=p.size();
  double phi=0;
//  #pragma omp parallel for reduction (+:phi)
  for (int i=0; i<j; i++)
    phi+=p[i].potential(p[j]);
//  #pragma omp parallel for reduction (+:phi)
  for (int i=j+1; i<ps; i++)
    phi+=p[i].potential(p[j]);
  return lB*p[j].charge*phi;
};
double interact::energy_dh(vector<particle> &p, int j) {
  if (p[j].charge==0) return 0;
  int ps=p.size();
  double u=0;
  for (int i=0; i<j; i++)
    u += energy_dh( p[i], p[j] );
  for (int i=j+1; i<ps; i++)
    u += energy_dh( p[i], p[j] );
  return lB*u;
};
/********************************
   GROUP WITH REST OF PARTICLES
 ********************************/
bool collision::overlap(vector<particle> &p, group &g) {
  int psize=p.size();
  if (g.beg!=-1) {
    for (int i=0; i<g.beg; i++)
      if (overlap(p, g, i)==true) return true;
    for (int i=g.end+1; i<psize; i++) {
      if (overlap(p, g, i)==true) return true;
    };
  };
  return false;
};
double interact::energy(vector<particle> &p, group &g) {
  int glen=g.end+1, psize=p.size();
  double u=0;
  if (g.beg!=-1)
    for (int i=g.beg; i<glen; i++) {
//      #pragma omp parallel for reduction (+:u)
      for (int j=0; j<g.beg; j++)
        u += energy(p[i],p[j]);
//      #pragma omp parallel for reduction (+:u)
      for (int j=glen; j<psize; j++)
        u += energy(p[i],p[j]);
    };
  return lB*u;
};

double interact::energy_dh(vector<particle> &p, group &g) { 
  int glen=g.end+1, psize=p.size(); 
  double u=0; 
  if (g.beg!=-1) 
    for (int i=g.beg; i<glen; i++) { 
      for (int j=0; j<g.beg; j++) 
        u += energy_dh(p[i],p[j]);                                                            
      for (int j=glen; j<psize; j++)                                                       
        u += energy_dh(p[i],p[j]);                                                            
    };                                                                                     
  return lB*u;                                                                             
};                                                                                         

/*******************************/
//group with j'th particle
bool collision::overlap(vector<particle> &p, group &g, int j) {
  if (g.beg==-1) return false;
  if (g.radius>0 && g.cm!=-1)
    if ( abs(p[g.cm].dist(p[j])) > g.radius+p[j].radius)
      return false;
  int len=g.end+1;

  if (g.isingroup(j)==false) {          //check if j is part of g
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
};

double interact::energy(vector<particle> &p, group &g, int j) {
  if (p[j].charge==0 || g.beg==-1) return 0;
  double u=0;
  int len=g.end+1;
  if (g.isingroup(j)==true) {   //avoid self-interaction...
//    #pragma omp parallel for reduction (+:u)
    for (int i=g.beg; i<j; i++)
      u+=energy(p[i],p[j]);
//    #pragma omp parallel for reduction (+:u)
    for (int i=j+1; i<len; i++)
      u+=energy(p[i],p[j]);
  } else                        //simple - j not in g
//    #pragma omp parallel for reduction (+:u)
    for (int i=g.beg; i<len; i++)
      u+=energy(p[i],p[j]);
  return lB*u;  
};

/********************************
 *
 * GROUP <-> ARBITRARY PARTICLE
 *
 ********************************/
bool collision::overlap(vector<particle> &p, group &g, particle &a) {
  int len=g.end+1;
  if (g.beg!=-1)
    for (int i=g.beg; i<len; i++)
      if (p[i].overlap(a)==true)
	return true;
  return false;
};

double interact::energy(vector<particle> &p, group &g, particle &a) {
  if (a.charge==0 || g.beg==-1) return 0;
  double phi=0;
  int len=g.end+1;
  if (g.beg!=-1)
//    #pragma omp parallel for reduction (+:phi)
    for (int i=g.beg; i<len; i++) 
      phi+=p[i].potential(a); 
  return lB*a.charge*phi;
};


/*********************
   SYSTEM ENERGY
 *********************/
double interact::energy(vector<particle> &p) {
  double u=0;
  double qq;
  unsigned int n = p.size();
  for (unsigned int i=0; i<n-1; i++)
    for (unsigned int j=i+1; j<n; j++)
      if (p[i].charge * p[j].charge != 0)
        u += energy(p[i], p[j]);
  return lB*u; // (kT)
};
double interact::energy_dh(vector<particle> &p) {
  double u=0;
  for (int i=0; i<p.size()-1; i++)
    for (int j=i+1; j<p.size(); j++)
      if (p[i].charge!=0 && p[j].charge!=0)
        u += energy_dh(p[i], p[j]);
  return lB*u;
};

/**********************
   BETWEEN TWO GROUPS
 **********************/
double interact::energy(vector<particle> &p, group &g1, group &g2) {
  if (g1.beg==-1 || g1.beg==-1)
    return 0.;

  unsigned int ilen=g1.end+1;                                 
  unsigned int jlen=g2.end+1;                                 
  if (g1.cm!=-1 && ilen>g1.cm) ilen=g1.cm;  //avoid ghost particles
  if (g2.cm!=-1 && jlen>g2.cm) jlen=g2.cm;  //such as "cm" and "dipv"
  
  double u=0;
  if (g1.vdw==true && g2.vdw==true && vdw!=0)
    for (int i=g1.beg; i<ilen; i++)
      for (int j=g2.beg; j<jlen; j++)
        u += energy_vdwchg(p[i],p[j]);
  else {
    for (int i=g1.beg; i<ilen; i++)
      for (int j=g2.beg; j<jlen; j++)
        u += energy(p[i],p[j]);
  };
  return lB*u;
};

double interact::energy_vdw(vector<particle> &p, group &g1, group &g2) {
  if (g1.beg==-1 || g2.beg==-1)
    return 0.;
  double u=0, r2;
  unsigned int ilen=g1.end+1;
  unsigned int jlen=g2.end+1;
  if (g1.cm!=-1 && ilen>g1.cm) ilen=g1.cm; // avoid ghost particles
  if (g2.cm!=-1 && jlen>g2.cm) jlen=g2.cm; // such as "cm" and "dipv" (see group.h)
  if (g1.vdw==true && g2.vdw==true && vdw!=0)
    for (int i=g1.beg; i<ilen; i++)
      for (int j=g2.beg; j<jlen; j++) {
        r2=p[i].sqdist(p[j]);
        u -= vdw / (r2*r2*r2); // in kT
      };
  return u; // (kT)
};

double interact::energy_dh(vector<particle> &p, group &g1, group &g2) {
  if (g1.beg==-1 || g2.beg==-1) return 0.;
  int ilen=g1.end+1, jlen=g2.end+1;
  double u=0;
  for (int i=g1.beg; i<ilen; i++)
    for (int j=g2.beg; j<jlen; j++)
      u += energy_dh(p[i],p[j]);
  return lB*u;
};

bool collision::overlap(vector<particle> &p, group &g1, group &g2) {
  if (g1.beg==-1 || g2.beg==-1) return false;
  if (g1.radius!=0 && g2.radius!=0)
    if (g1.cm!=-1 && g2.cm!=-1)
      if (abs(p[g1.cm].dist(p[g2.cm])) > g1.radius+g2.radius)
        return false;
  int ilen=g1.end+1, jlen=g2.end+1;
  for (int i=g1.beg; i<ilen; i++) {
    for (int j=g2.beg; j<jlen; j++) {
      if ( p[i].overlap(p[j])==true )
        return true;
    };
  };
  return false;
};

//between first in list "groups" and the remaining of "groups"
double interact::energy(vector<particle> &p, vector<group> &g, int groups,...) {
  int j,i=groups;
  double u=0;
  va_list ap;
  va_start(ap, groups);
  for (;;) {
    j=va_arg(ap, int);
    if (j==-1)
      break;
    u+=energy( p, g[i], g[j] ); //already in kT
  };
  return u; //do not multiply with lB!
};

//non-electrostatic energy of chain particle i with rest of the chain
double interact::chain(vector<particle> &p, group &g, int i) {
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
};

//graft energy of the chain, g.
double interact::graft(vector<particle> &p, group &g) {
  if (g.graftpoint!=-1)
    return quadratic( p[g.graftpoint], p[g.beg]);
  else
    return 0;
};

/*!
 * ...between the two dipoles a and b, separated by the
 * distance r.
 * \f$ \beta u(r) = l_B \frac{a_x b_x + a_y b_y - 2a_z b_z  }{r^3}\f$
 */
double interact::dipdip(point &a, point &b, double r) {
  return lB*( a.x*b.x + a.y*b.y - 2*a.z*b.z )/(r*r*r);
};
double interact::iondip(point &a, double q, double r) {
  return -lB*q*a.z/(r*r);
};

//Metropolis test. du must be in units of kT!
bool interact::metropolis(double du) {
  if (du > 0)
    if ( slmp.random_one()>exp(-du) )
      return false;
  return true;
};

// Total electrostatic potential in a point
double interact::potential(vector<particle> &p, point &a) {
  double u=0;
  int ps=p.size();  
  for (int i=0; i<ps; i++)
    u+=p[i].charge * p[i].invdist(a);
  return lB*u;
};

// Internal (NON)-electrostatic energy in group
double interact::internal(vector<particle> &p, group &g) {
  int glen=g.end+1, psize=p.size();
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
        u+=energy(p[i],p[j]);
  }
  return lB*u;
};

// Internal energy of particles specified in an vector
double interact::internal(vector<particle> &p, vector<short int> &g) {
  double u=0;
  int size=g.size();
  for (int i=0; i<size-1; i++)
    for (int j=i+1; j<size; j++)
      u+=energy(p[g[i]], p[g[j]]);
  return lB*u;
};
