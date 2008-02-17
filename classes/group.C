#include "group.h"

//---------------- GROUP -----------------------
group::group(int i) {
  set(-1,-1);     //is empty
  title="GROUP";
  name="(noname)";//is nameless
}

void group::set(short int first, short int last) {
  beg=first;
  end=last;
}

short int group::size() { return (beg==-1) ? 0 : end-beg+1; }

/*! \brief Calculates total charge
 *  \return Charge number. The charge is also stored in
 *          cm.charge.
 */
double group::charge(vector<particle> &p) {
  double z=0;
  for (short i=beg; i<=end; i++)
    z+=p[i].charge;
  cm.charge=z;
  return z;
}

void group::operator+=(group g) {
  if (g.beg==-1 && g.end==-1)   // if added group is empty
    return;
  if (beg==-1 && end==-1)       // if this is empty
    beg=g.beg;
  end=g.end;
}

group group::operator+(group g) {
  group o;
  if (beg<g.beg) {
    o.beg=beg;
    o.end=g.end;
  } else {
    o.beg=g.beg;
    o.end=end;
  };
  o.name = this->name + " + " + g.name;
  if (o.size()!=this->size()+g.size())
    cout << "# Warning: Added groups are not continous!\n";
  return o;
}

short int group::random() { return beg + rand() % size(); }

string group::info() {
  ostringstream o;
  o << endl
    << "# " << title << ": " << name << endl
    << "#   Range                  = [" << beg << "," << end << "]" << " " << size() << endl
    << "#   Mass center            = " << cm << endl;
  return o.str();
}

bool group::find(unsigned int i) { return (i>=beg && i<=end) ? true : false; }

/*!
 * Recalcute mass center of the particles
 * in a group. The result is stored in the
 * group CM placeholder as well as returned
 * by the function.
 * \warning Doesn't consider periodic boundaries!
 * \note Use masscenter(container &) instead.
 */
point group::masscenter(vector<particle> &p) {
  cm.clear();
  double sum=0;
  for (short i=beg; i<=end; i++) {
    cm += p[i] * p[i].mw;
    sum += p[i].mw;
  }
  cm=cm*(1./sum);
  cm_trial = cm;
  return cm;
}

/*!
 * This mass-center routine obeys periodic boundaries (if only)
 * by first centering the first particle in (0,0,0), calc. CM
 * and move it back.
 *
 * \author Mikael Lund
 * \date Dec. 2007, Prague
 */
point group::masscenter(container &con) {
  double sum=0;
  cm.clear();
  point t, o = con.p[beg]; // set origo to first particle
  for (unsigned short i=beg; i<=end; i++) {
    t = con.p[i]-o;        // translate to origo
    con.boundary(t);       // periodic boundary (if any)
    cm += t * con.p[i].mw;
    sum += con.p[i].mw; 
  }
  cm=cm*(1./sum) + o;
  con.boundary(cm);
  cm_trial = cm;
  return cm;
}


void group::undo(particles &par) {
  cm_trial = cm;
  for (short i=beg; i<=end; i++) {
    par.trial[i].x = par.p[i].x;
    par.trial[i].y = par.p[i].y;
    par.trial[i].z = par.p[i].z;
  }
}
void group::accept(particles &par) {
  cm = cm_trial;
  unsigned short i,ilen=end+1;
  for (i=beg; i<ilen; i++) {
    par.p[i].x = par.trial[i].x;
    par.p[i].y = par.trial[i].y;
    par.p[i].z = par.trial[i].z;
  }
}
void group::add(container &par, vector<particle> v, bool collision) {
  beg=par.p.size();
  for (unsigned short i=0; i<v.size(); i++) {
    par.boundary(v[i]);
    par.push_back(v[i]);
    end=beg+i;
  }
  masscenter(par);                // calc. mass center
  
  // test for overlap w. other particles
  if (collision==true) {
    move(par, -cm);               // translate to origo (0,0,0)
    accept(par); 
    point a;
    while (overlap(par)==true) {  // No overlap allowed
      par.randompos(a);           // find a random, vacent 
      move(par, a);               // position.
      accept(par);
    }
  }
}
void group::add(container &con, particle::type id, short n) {
  particle a=con.get(id);
  if (beg<0)
    beg=con.p.size();
  for (unsigned short i=0; i<n; i++) {
    do con.randompos(a);
    while (con.overlap(a)==true); 
    end = con.push_back(a)-1;
  }
}

unsigned short group::displace(container &c, double dp) {
  unsigned short i=random();
  c.trial[i].x = c.p[i].x + dp*slp.random_half();
  c.trial[i].y = c.p[i].y + dp*slp.random_half();
  c.trial[i].z = c.p[i].z + dp*slp.random_half();
  c.boundary(c.trial[i]);
  return i;
}

/*!
 * \param par Container class
 * \param c ...by adding this vector to all particles
 * \param k Keyword to specify if the move should be automatically accepted (default no).
 */
void group::move(container &par, point c) {
  for (unsigned short i=beg; i<=end; i++) {
    par.trial[i].x = par.p[i].x + c.x;
    par.trial[i].y = par.p[i].y + c.y;
    par.trial[i].z = par.p[i].z + c.z;
    par.boundary(par.trial[i]);
  }
  cm_trial.x = cm.x + c.x;
  cm_trial.y = cm.y + c.y;
  cm_trial.z = cm.z + c.z;
  par.boundary(cm_trial);
}

bool group::overlap(container &c) {
  double s;
  unsigned short i,j;
  for (i=beg; i<=end; i++) {
    if (c.collision(c.p[i])==true)  // check cell collision
      return true;
    for (j=0; j<beg; j++) {
      s=c.p[i].radius+c.p[j].radius;
      if ( c.sqdist(c.p[i],c.p[j]) < s*s )
        return true;
    }
    for (j=end+1; j<c.p.size(); j++) {
      s=c.p[i].radius+c.p[j].radius;
      if ( c.sqdist(c.p[i],c.p[j]) < s*s )
        return true;
    }
  }
  return false;
}
short int group::count(vector<particle> &p, particle::type id) {
  short int i,n=0;
  for (i=0; i<p.size(); i++)
    if (p[i].id==id)
      n++;
  return n;
}
unsigned short group::nummolecules() { return size(); }

/*****************************
 *
 *          S A L T
 *
 *****************************/
salt::salt(particle::type cat, particle::type an)
{
  title="MOBILE SALT PARTICLES";
  cation=cat;
  anion=cat;
}
string salt::info(container &con) {
  ostringstream o;
  float c=1./6.022e23/1e-27;
  short nan=count(con.p, anion),
        ncat=count(con.p, cation);
  o << group::info();
  o << "#   Cation (type,z,N,conc) = "
    << con.d[cation].name << ", " << con.d[cation].p.charge << ", "
    << ncat << ", " << ncat/con.getvolume()*c << endl
    << "#   Anion (type,z,N,conc)  = "
    << con.d[anion].name << ", " << con.d[anion].p.charge << ", "
    << nan << ", " << nan/con.getvolume()*c << endl;
  return o.str();
}
void salt::isobaricmove(container &con, double newlen) {
  for (unsigned short i=beg; i<=end; i++)
    con.scale(con.trial[i], newlen);
}

/*!
 * This searches the inputfile object for the
 * keywords "nion#" and "tion#" and if found
 * tries to insert the salt particles at random
 * positions.
 */
void salt::add(container &con, inputfile &in) {
  short n=1, npart;
  particle::type id;
  while (n<3) {
    ostringstream nion, tion;
    nion << "nion" << n;
    tion << "tion" << n++;
    npart = in.getint(nion.str(), 0);
    id = con.id(in.getstr(tion.str()));
    if (npart!=0)
      group::add(con, id, npart ); // add particles
    if (con.d[id].p.charge>0)
      cation=id;
    if (con.d[id].p.charge<0)
      anion=id;
  }
}

//--------------- MACROMOLECULE ---------------
macromolecule::macromolecule() { title="MACROMOLECULE"; }
string macromolecule::info() {
  ostringstream o;
  o << group::info();
  if (Q.cnt>0)
    o << "#   Charge <Z>             = " << Q.avg() << endl
      << "#   Charge Sqr <Z^2>       = " << Q2.avg() << endl
      << "#   Capacitance            = " << Q2.avg()-pow(Q.avg(),2) << endl;
  if (dip.cnt>0)
    o << "#   Dipole moment          = " << dip.avg() << " " << mu << endl;
  return o.str();
}

unsigned short macromolecule::nummolecules() { return 1; }

void macromolecule::operator=(group g) {
  beg=g.beg;
  end=g.end;
  name=g.name;
  cm=g.cm;
  cm_trial=g.cm_trial;
}
/*!
 * Recalulates the dipole moment of a
 * group and stores the result in the
 * group dipole moment placeholder.
 * Origo is (0,0,0).
 *
 * \todo Implement PBC!
 */
double macromolecule::dipole(vector<particle> &p)
{
  mu.clear();
  double a;
  for (unsigned short i=beg; i<=end; ++i) {
    mu.x += (p[i].x-cm.x) * p[i].charge;
    mu.y += (p[i].y-cm.y) * p[i].charge;
    mu.z += (p[i].z-cm.z) * p[i].charge;
  }
  a=mu.len();
  dip+=a;
  return a;
}
double macromolecule::radius(vector<particle> &p) {
  double r,max=0;
  masscenter(p);
  for (short i=beg; i<=end; i++) {
    r=p[i].dist(cm) + p[i].radius;
    if (r>max)
      max=r;
  }
  cm.radius = max;
  return max;
} 
double macromolecule::charge(vector<particle> &p) {
  double z=group::charge(p);
  Q+=z;
  Q2+=z*z;
  return z;
}

void macromolecule::center(container &con) {
  move(con, -cm);
  accept(con);
};

void macromolecule::zmove(container &par, double dz) {
  cm_trial.z = cm.z + dz;
  par.boundary(cm_trial);
  for (short int i=beg; i<=end; i++) {
    par.trial[i].z = par.p[i].z + dz;
    par.boundary(par.trial[i]);
  };
}
void macromolecule::rotate(container &par, double drot) {
  point u;
  double r=2;
  while (r > 1.) { //Generate a random unit vector
    u.x=slp.random_one();
    u.y=slp.random_one();
    u.z=slp.random_one();
    r=sqrt(u.x*u.x+u.y*u.y+u.z*u.z);
  }
  u.x=u.x/r;
  u.y=u.y/r; 
  u.z=u.z/r;
  rotate(par, u, drot*slp.random_half());
}

/*!
 * \param g Group to rotate
 * \param u Vector to rotate around
 * \param angle ..by this many degrees (rad)
 * \param k Keyword to specify automatic acceptance
 * \toto This could be done more elegant using quaternion parameters (Frenkel p49)
 */
void macromolecule::rotate(container &par, point u, double angle) {
  point b;
  double cosang, sinang, eb;
  double e1mcox, e1mcoy, e1mcoz;
  
  cosang=cos(angle);
  sinang=sin(angle);
  e1mcox=(1.-cosang)*u.x;
  e1mcoy=(1.-cosang)*u.y;
  e1mcoz=(1.-cosang)*u.z;
  for (short i=beg; i<=end; i++) {
    b.x=par.p[i].x - cm.x;              // translate to origo...
    b.y=par.p[i].y - cm.y;
    b.z=par.p[i].z - cm.z;
    par.boundary(b);                    // Apply boundary conditions
    eb=u.x*b.x + u.y*b.y + u.z*b.z;
    par.trial[i].x=e1mcox*eb+cosang*b.x+sinang*(u.y*b.z - u.z*b.y) + cm.x;
    par.trial[i].y=e1mcoy*eb+cosang*b.y+sinang*(u.z*b.x - u.x*b.z) + cm.y;
    par.trial[i].z=e1mcoz*eb+cosang*b.z+sinang*(u.x*b.y - u.y*b.x) + cm.z;
    par.boundary(par.trial[i]);
  }
}

void macromolecule::add(container &con, inputfile &in ) { }

void macromolecule::isobaricmove(container &con, double newlen) {
  double oldvol=con.getvolume(); // store original volume
  con.scale(cm_trial, newlen);   // scale mass center
  con.setvolume(pow(newlen,3));  // the boundary function needs the trial volume...
  move(con, cm_trial-cm );       // apply periodic boundaries
  con.setvolume(oldvol);         // restore original volume
}

//--------------- CHAIN -----------------
chain::chain() { graftpoint=-1; }

double chain::monomerenergy(vector<particle> &p, short i) {
  double u=0;
  //the first ?
  if (i==beg) {
    u+=quadratic( p[i], p[i+1] );
    if (graftpoint>-1)
      u+=quadratic( p[i], p[graftpoint] );
    return u;
  }

  //the last ?
  if (i==end)
    return quadratic( p[i], p[i-1] );

  //otherwise...
  return quadratic(p[i], p[i+1]) + quadratic(p[i], p[i-1]);
}

double chain::internalenergy(vector<particle> &p) {
  double u=0;
  for (short i=beg; i<end; i++)
    u+=quadratic( p[i], p[i+1]);
  if (graftpoint>-1)
    u+=quadratic( p[beg], p[graftpoint] );
  return u;
}

//-------------------- PLANE -----------------
//! By default the plane lies in the XY plane at
//! z=0. This plane can be put in XZ and YZ by changing
//! the plane vector. To shift the plane, modify the
//! offset vector.
planarsurface::planarsurface() {
  plane.x=1;
  plane.y=1;
  plane.z=0;
}
void planarsurface::project(point &p) { p=p*plane + offset; }

//--------------------- PHOSPHOMEMBRANE -------------------
double zwittermembrane::selfenergy(particles &par) {
  double u=0;
  for (short i=beg; i<end; i=i+2)
    u+=pairpot(par.p[i], par.p[i+1]);
  return u;
}
short zwittermembrane::mate(short i) { return (i%2==0) ? i+1 : i-1; }

//! Even particles can move only in the plane while odd
//! particles can move in any direction. The latter
//! corresponds to flexible N(CH3)+ groups attached to
//! phosphate groups. 
unsigned short zwittermembrane::displace(container &c, double d) {
  unsigned short i=group::displace(c,d);
  if (i%2==0)
    project(c.trial[i]); // keep even particles in the plane
  return i;
}       
