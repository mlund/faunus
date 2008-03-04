#include "analysis.h"

//! \param u_init Initial total system energy
systemenergy::systemenergy(double u_init) {
  u0=u_init;
  sum=u0;
  cur=u0;
}
void systemenergy::update(double energy) {
  cur=energy;
  uavg+=cur;
  u2avg+=cur*cur;
}
void systemenergy::operator+=(double du) { sum+=du; }

//---------------- ANGULAR CORRELATIONS ---------------------------
angularcorr::angularcorr() {}
void angularcorr::update(macromolecule &g1, macromolecule &g2, distributions &d) {
  double l1=g1.mu.len(), l2=g2.mu.len(), z;
  if (l1>0 && l2>0 )  {
    m1=g1.mu*(1./l1);
    m2=g1.mu*(1./l2); // unit vectors at origo
    z=g1.cm.dist(g2.cm);
    d.add("dip_z*z", z, m1.z*m2.z);
    d.add("dip_x*x", z, m1.x*m2.x);
    d.add("dip_y*y", z, m1.y*m2.y);
  }
}

//---------------- DISTRIBUTIONS ----------------------------------
/*!
 * By default the maximum and minimum x values are dynamically
 * resized to match any given dataset (for x>=0).
 *
 * \note If x<0, "min" should be specified
 * \param deltax x value resolution
 * \param min Minimum x value (optional, auto-resizes)
 * \param max Maximum x value (optional, auto-resizes)
 */
distributions::distributions(float deltax, float min, float max) {
  dx=deltax;
  xmax_set=max;
  xmin_set=min;
  xmax=-1e7;
  xmin=-xmax;
}

/*! 
 * Adds a value to the distribution. The y value is automatically
 * averaged. If "name" is not found in the current distribution,
 * it will be created and the data point will be added.
 * 
 * \param name Name of distribution to add to.
 * \param x x-value
 * \param y value (add to average)
 */
bool distributions::add( string name, float x, float y )
{
  unsigned short i = find(name);
  d[i](x)+=y; // add to average
  if (x>xmax) // track max x value
    xmax=x;
  if (x<xmin) // track min x value
    xmin=x;
  return true;
}

unsigned short distributions::find( string name ) {
  unsigned short i;
  for (i=0; i<s.size(); i++) // search name vector
    if (s[i]==name)
      return i;
  s.push_back(name); // if not found, add it!
  d.resize( d.size()+1,
      xytable<float, average<float> >(dx,xmin_set,xmax_set));
  return s.size()-1;
}

string distributions::info() {
  unsigned char i;
  ostringstream o;
  o << "# DISTRIBUTION FUNCTIONS:\n"
    << "# 1 = distance" << endl;
  for (i=0; i<s.size(); i++)
    o << "# " << i+2 << " = " << s[i] << endl;
  for (float x=xmin; x<=xmax; x+=dx) {
    o << x;
    for (i=0; i<d.size(); i++)
      o << " " << d[i](x).avg();
    o << endl;
  }
  return o.str();
}

bool distributions::write(string name)
{
  return fio.writefile(name, info() );
}
