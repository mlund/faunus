#include "particles.h"

int particles::push_back(particle &par) {
  p.push_back(par);
  trial.push_back(par);
  return p.size();
};

group particles::append(vector<particle> v) {
  group g;
  g.beg=p.size();
  for (unsigned int i=0; i<v.size(); i++) {
    push_back(v[i]);
    g.end=g.beg+i;
  };
  return g;
};

double particles::charge() {
  double z=0;
  for (int i=0; i<p.size(); i++)
    z+=p[i].charge;
  return z;
}

/*!\param origo Center of the spherical region
 * \param r Radius of the spherical region
 */
double particles::charge(point &origo, double r) {
  double q=0,r2=r*r;
  for (int i=0; i<p.size(); i++)
    if (p[i].sqdist(origo) <= r2)
      q += p[i].charge;
  return q;
};

double particles::radius(group &g, point &origo, keys k) {
  double d, max=0,min=10000.;
  average<double> r;
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

void particles::zmove(group &g, double dz, keys k) {
  for (int i=g.beg; i<=g.end; i++)
    trial[i].z = p[i].z + dz;
  g.cm_trial.z = g.cm.z + dz;
  if (k==AUTOACCEPT)
    accept(g);
}

/*!
 * \param g Group to translate...
 * \param c ...by adding this vector to all particles
 * \param k Keyword to specify if the move should be automatically accepted (default no).
 */
void particles::move(group &g, point c, keys k) {
  for (int i=g.beg; i<=g.end; i++) {
    trial[i].x = p[i].x + c.x;
    trial[i].y = p[i].y + c.y;
    trial[i].z = p[i].z + c.z;
  };
  g.cm_trial = g.cm + c;

  if (k==AUTOACCEPT) accept(g);
};
void particles::displace(int i, double dp) {
  trial[i].x = p[i].x + dp*slp.random_half();
  trial[i].y = p[i].y + dp*slp.random_half();
  trial[i].z = p[i].z + dp*slp.random_half();
}

void particles::rotate(group &g, double drot, keys k) {
  point u;
  double r=2;
  while (r > 1.) { //random unit vector
    u.x=slp.random_one();
    u.y=slp.random_one();
    u.z=slp.random_one();
    r=sqrt(u.x*u.x+u.y*u.y+u.z*u.z);
  };
  u.x=u.x/r;
  u.y=u.y/r; 
  u.z=u.z/r;
  rotate(g, u, drot*slp.random_half(), k );
}

/*!
 * \param g Group to rotate
 * \param u Vector to rotate around
 * \param angle ..by this many degrees (rad)
 * \param k Keyword to specify automatic acceptance
 */
void particles::rotate(group &g, point u, double angle, keys k) {
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
    bz=p[i].z - g.cm.z; //p[g.cm].z;
    eb=u.x*bx + u.y*by + u.z*bz;
    trial[i].x=e1mcox*eb+cosang*bx+sinang*(u.y*bz - u.z*by);
    trial[i].y=e1mcoy*eb+cosang*by+sinang*(u.z*bx - u.x*bz);
    trial[i].z=e1mcoz*eb+cosang*bz+sinang*(u.x*by - u.y*bx) + g.cm.z; //p[g.cm].z;
  };
  if (k==AUTOACCEPT) accept(g);
}

void particles::undo(group &g, keys k) {
  int ilen=g.end+1;
  g.cm_trial = g.cm;
  if (k==ALL)
    for (int i=g.beg; i<ilen; i++)
      trial[i]=p[i];
  else
    for (int i=g.beg; i<=g.end; i++) {
      trial[i].x = p[i].x;
      trial[i].y = p[i].y;
      trial[i].z = p[i].z;
    };
}

void particles::accept(group &g, keys k) {
  int ilen=g.end+1;
  g.cm = g.cm_trial;
  if (k==ALL)
    for (int i=g.beg; i<ilen; i++)
      p[i]=trial[i];
  else
    for (int i=g.beg; i<ilen; i++) {
      p[i].x = trial[i].x;
      p[i].y = trial[i].y;
      p[i].z = trial[i].z;
    };
}

