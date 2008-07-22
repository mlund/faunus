#include "space.h"

// VECTOR HANDLING
int space::push_back(particle &par) {
  p.push_back(par);
  trial.push_back(par);
  return p.size();
};

int space::trialpush_back(particle &par) {
  trial.push_back(par);
  return trial.size();
};

//ADD A PARTICLE VECTOR
group space::append(vector<particle> v) {
  group g;
  g.beg=p.size();
  for (int i=0; i<v.size(); i++) {
    push_back(v[i]);
    g.end=g.beg+i;
  };
  return g;
};


// PROPERTIES
double space::charge(group &g) {
  double z=0;
  for (int i=g.beg; i<g.end+1; i++)
    z+=p[i].charge;
  g.charge=z;
  return z;
};
double space::charge() {
  double z=0;
  for (int i=0; i<p.size(); i++)
    z+=p[i].charge;
  return z;
}
double space::charge(point &origo, double r) {
  double q=0,rsq=r*r;
  for (int i=0; i<p.size(); i++)
    if (p[i].sqdist(origo) <= rsq)
      q += p[i].charge;
  return q;
};

double space::radius(group &g, point &origo, keys k) {
  double d, max=0,min=10000.;
  average r;
  switch (k) {
    case MAX :
      for (int i=g.beg; i<=g.end; i++) {
        d=p[i].dist(origo) + p[i].radius;
        if (d>max)
          max=d;
      };
      return max;
      break;
    case AVERAGE :
      for (int i=g.beg; i<=g.end; i++)
        r += p[i].dist(origo) + p[i].radius;
      return r.avg();
      break;
    case INTERIOR :
      for (int i=g.beg; i<=g.end; i++)
        if (p[i].charge!=0) {
          d=p[i].dist(origo);
          if (d<min) min=d;
        };
      return min;
      break;
  };
};


//calculate center-of-mass
point space::center_of_mass(group &g) {
  double sum=0;
  point c;
  for (int i=g.beg; i<=g.end; i++) {
    c += p[i]*p[i].mw;
    sum += p[i].mw;
  };
  c = c * (1./sum);
  return c;
};

void space::addmasscenter(group &g) {
  particle cm;
  cm.mw=1e-6;
  cm=center_of_mass(g);
  push_back(cm);        // add to end of vectors (trial and p)
  g.cm=p.size()-1;      // store its position
  if (g.cm==g.end+1)
    g++;
  else
    cout << "# Warning: Mass center outside group range!" << endl;
};

//Sort so that charged species are in front
group space::sort(group &g) {
  group out = g;
  int j=0;
  vector<particle> s(0);
  for (int i=g.beg; i<=g.end; i++) //find charges
    if (p[i].charge != 0) s.push_back(p[i]);

  out.end = out.beg + s.size() - 1; //end of charges
  
  for (int i=g.beg; i<=g.end; i++) //find neutral
    if (p[i].charge == 0) s.push_back(p[i]);
  for (int i=g.beg; i<=g.end; i++) //copy into particle vector
    p[i]=s[j++];

  return out;
};

/***************************************
 *                                     *
 * TRANSLATIONAL AND ROTATIONAL MOVES  *
 *                                     *
 ***************************************/
void space::zmove(group &g, double dz, keys k) {
  for (int i=g.beg; i<=g.end; i++)
    trial[i].z = p[i].z + dz;
  if (k==AUTOACCEPT) accept(g);
};
void space::move(group &g, point c, keys k) {
  for (int i=g.beg; i<=g.end; i++) {
    trial[i].x = p[i].x + c.x;
    trial[i].y = p[i].y + c.y;
    trial[i].z = p[i].z + c.z;
  };
  if (k==AUTOACCEPT) accept(g);
};
void space::displace(int i, double dp) {
  trial[i].x = p[i].x + dp*random_half();
  trial[i].y = p[i].y + dp*random_half();
  trial[i].z = p[i].z + dp*random_half();
};

void space::rotate(group &g, double drot, keys k) {
  point u;
  double r=2;
  while (r > 1.) { //random unit vector
    u.x=random_one();
    u.y=random_one();
    u.z=random_one();
    r=sqrt(u.x*u.x+u.y*u.y+u.z*u.z);
  };
  u.x=u.x/r;
  u.y=u.y/r; 
  u.z=u.z/r;
  rotate(g, u, drot*random_half(), k );
};

void space::rotate(group &g, point u, double angle, keys k) {
  double bx, by, bz;
  double cosang, sinang, eb;
  double e1mcox, e1mcoy, e1mcoz;
  
  cosang=cos(angle);
  sinang=sin(angle);
  e1mcox=(1.-cosang)*u.x;
  e1mcoy=(1.-cosang)*u.y;
  e1mcoz=(1.-cosang)*u.z;

  for (int i=g.beg; i<=g.end; i++) {
    bx=p[i].x;
    by=p[i].y;
    bz=p[i].z - p[g.cm].z;
    eb=u.x*bx + u.y*by + u.z*bz;
    trial[i].x=e1mcox*eb+cosang*bx+sinang*(u.y*bz - u.z*by) ;
    trial[i].y=e1mcoy*eb+cosang*by+sinang*(u.z*bx - u.x*bz) ;
    trial[i].z=e1mcoz*eb+cosang*bz+sinang*(u.x*by - u.y*bx) + p[g.cm].z;
  };
  if (k==AUTOACCEPT) accept(g);
};


// UNDO FUNCTIONS
void space::undo(group &g, keys k) {
  int ilen=g.end+1;
  if (k==ALL)
    for (int i=g.beg; i<ilen; i++)
      trial[i]=p[i];
  else
    for (int i=g.beg; i<=g.end; i++) {
      trial[i].x = p[i].x;
      trial[i].y = p[i].y;
      trial[i].z = p[i].z;
    };
};

void space::accept(group &g, keys k) {
  int ilen=g.end+1;
  if (k==ALL)
    for (int i=g.beg; i<ilen; i++)
      p[i]=trial[i];
  else
    for (int i=g.beg; i<ilen; i++) {
      p[i].x = trial[i].x;
      p[i].y = trial[i].y;
      p[i].z = trial[i].z;
    };
};

// DATA HANDLING AND MANIPULATION
group space::insert_chain(chain::chain &c) {
  group g;
  collision col;

  g.graftpoint=c.graftpoint;
  g.beg=p.size();
  g.chain=true;

  //put near graftpoint (if any).
  if (c.graftpoint > -1) {
    c.v[0].x = p[c.graftpoint].x;
    c.v[0].y = p[c.graftpoint].y;
    c.v[0].z = p[c.graftpoint].z;
  } else {
    c.v[0].x = 40.*random_half();
    c.v[0].y = 40.*random_half();
    c.v[0].z = 0;
  };

  //insert chain...
  for (int i=0; i<c.v.size(); i++) {

    //place monomer over former
    if (i!=0) {
      c.v[i].x=p[g.end-1].x;
      c.v[i].y=p[g.end-1].y;
      c.v[i].z=p[g.end-1].z;
    };

    //add it
    g.end=push_back( c.v[i] )-1;

    //keep on if overlap
    while (col.overlap(p, g.end)==true) {
      p[g.end]=c.v[i];
      displace( g.end, 10. ); // <- notice displacement parameter.
      p[g.end]=trial[g.end];
    };    
    g.end++;
  };
  g.end=p.size()-1;
  //cout << "Inserted!\n";
  return g;
};

/*! \fn group space::insert_salt( int n, double charge, double radius, cell &c)
 *  \brief Adds mobile ions
 *  \param n Number of particles
 *  \param charge Valency of particles
 *  \param radius Radius of particles
 */
group space::insert_salt(int n, double charge, double radius, cell &c) {
  collision col;
  group g;
  if (n==0) return g;

  g.beg=p.size();

  particle a;
  a.radius=radius;
  a.charge=charge;
  
  for (int i=0; i<n; i++) {
    g.end=push_back(a)-1;
  //  c.randomPos(p[g.end]);      //Unchecked change
    while (col.overlap(p, g.end)==true) {
      c.randomPos(p[g.end]);
      trial[g.end] = p[g.end];
    };
  };
  return g;
};

bool space::save(string file) {
  ofstream f( file.c_str() );
  if (f) {
    f.precision(30);
    for (int i=0; i<p.size(); i++) {
      f << p[i].x << " "
        << p[i].y << " "
        << p[i].z << " "
        << p[i].charge << endl;
    };
    f.close();
    return true;
  };
  cout << "# Error writing state!!\n";
  return false;
};

// Read saved coordinates from disk
bool space::load(string file) {
  ifstream f( file.c_str() );
  if (f) {
    f.precision(30);
    for (int i=0; i<p.size(); i++)
      f >> p[i].x >> p[i].y >> p[i].z >> p[i].charge;
    f.close();
    trial = p;
    cout << "# Previously saved configuration loaded.\n";
    return true;
  };
  cout << "# Initial configuration not present.\n";
  return false;
};


//test if 'p' and 'trial' are idential
// (as they should right after accepted or rejected)
// (true if test is OK)
bool space::safetytest_vector() {
  bool rc=false;
  if (p.size()==trial.size()) {
    for (int i=0; i<p.size(); i++) {
      if (p[i].x!=trial[i].x ||
	  p[i].y!=trial[i].y ||
	  p[i].z!=trial[i].z || 
	  p[i].charge!=trial[i].charge)
      {
        cout << "#p[" << i <<"] miss matched"<< endl; 
        cout << "#p("<< p[i].x <<","<<p[i].y <<","<<p[i].z<<") q="<< p[i].charge<<endl;
        cout << "#trial("<< trial[i].x <<","<<trial[i].y <<","<<trial[i].z<<") q="<< trial[i].charge<<endl;

	rc=false;
	break;
      } else
        rc=true;
    };

  } else { cout << "Size dont match"<<endl;
  }
  if (rc==false)
    cout << "# Fatal error: Particle vectors corrupted!!\n";
  return rc;
};


/*------------------------------
  Dipole moment unit vector.

  g = group
  p = particle vector
  ------------------------------*/
point space::calcdipole(group &g)
{
  point mu;
  for (int i=g.beg; i<=g.end; i++) {
    //mu += ( p[i]-p[g.cm] ) * p[i].charge;
    mu.x += (p[i].x-p[g.cm].x) * p[i].charge;
    mu.y += (p[i].y-p[g.cm].y) * p[i].charge;
    mu.z += (p[i].z-p[g.cm].z) * p[i].charge;
  };
  
  g.dip = mu.len();

  if (g.dip!=0) {
    //mu = mu * (1./g.dip) + p[g.cm];
    mu.x = mu.x / g.dip;
    mu.y = mu.y / g.dip;
    mu.z = mu.z / g.dip;
  };
  mu.x += p[g.cm].x;
  mu.y += p[g.cm].y;
  mu.z += p[g.cm].z;
    
  p[g.dipv] = mu;
  trial[g.dipv] = mu;
  return (mu-p[g.cm])*g.dip ;
};

bool space::adddipole(group &g)
{
  particle mu;
  push_back(mu);         // add dipole to end of particle vectorS

  g.dipv=p.size()-1;     // set pointer to dipole particle
  if (g.dipv==g.end+1) {   // check range!
    g++;                 // expand range
    calcdipole(g);       // calculate dipole moment
    return true;
  };
  cout << "# warning: Mass center outside group range!" << endl;
  return false;
};

int space::cntHydrophobic(group &g)
{
  int cnt=0;
  for (int i=g.beg; i<=g.end; i++)
    if (p[i].hydr==true)
      cnt++;
  return cnt;
};

double space::rho(point &p0, group &g, double r, double Z, double delta) {
  double d,vol,n=0;
  for (int i=g.beg; i<=g.end; i++)
    if (p[i].charge==Z) {
      d=p0.dist(p[i]);
      if (d>=r && d<r+delta) n+=1.;
    };
  vol=4./3.*3.14 * ( pow(r+delta,3)-r*r*r );
  return n / vol;
};

int space::count(group &g, double Z) {
  int j=0;
  for (int i=g.beg; i<=g.end; i++)
    if (p[i].charge==Z) j++;
  return j;
};

int space::trialcount(group &g, double Z) {
  int j=0;
  for (int i=g.beg; i<=g.end; i++)
    if (trial[i].charge==Z) j++;
  return j;
};

void space::compoChangeIn(vector<particle> &pair, group &salt, keys k) {
  for (int i=0; i<2; i++ ) {
    trial.push_back(pair[i]);
  };
  salt.end=salt.end+2;
};

vector<int> space::compoChangeOut(group &salt, keys k) {
  vector<int> pair(2);
  int S=salt.size();
  int n, m;
  m=salt.beg;
  n=salt.end;
  while (trial[n].charge != 1)
    n=m+S*random_one();
  pair[0]=n;
  while (trial[n].charge != -1)
    n=m+S*random_one();
  pair[1]=n;
  salt.end=salt.end-1;
//  if (k==AUTOACCEPT)
//    accept(g);
  return pair;
};

