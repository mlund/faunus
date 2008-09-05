#include "analysis.h"

namespace Faunus {

//! \param u_init Initial total system energy
systemenergy::systemenergy(double u_init) {
  u0=u_init;
  sum=u0;
  cur=u0;
  confu.clear();
}
void systemenergy::update(double energy) {
  cur=energy;
  uavg+=cur;
  u2avg+=cur*cur;
}
void systemenergy::track() {
  confu.push_back(sum);
}
void systemenergy::write() {
  if (confu.size())
    fio.writefile("systemenergy.dat", confuout() );
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

//aggregation
aggregation::aggregation(container &C, vector<macromolecule> &G, double s) {
  con=&C;
  g.clear();
  for (int iter=0;iter<G.size();iter++)
    g.push_back(&(G[iter]));
  //g=&G;
  CNT=0;
  dist.clear();
  dist.resize(g.size(), 0);
  RG2.clear();
  average<double> start;
  for (int i=0;i<g.size();i++)
  RG2.push_back(start);
  sep=s;  
}
void aggregation::count() {
  ++CNT;
  for (int k=0; k<g.size();k++) {  //start looking for aggregates around k
  //prepare the vectors
    agg.clear();
    unagg.clear();
    for (int iter=0;iter<g.size();iter++)
      if (iter!=k)
        unagg.push_back(g[iter]);
    agg.push_back(g[k]);
  //start looking
    for (int i=0; i<agg.size(); i++){
      for (int j=0; j<unagg.size(); j++)
        if (coll.overlap(con->p, *agg[i], *unagg[j], sep) ) // Find the aggregated
          agg.push_back(unagg[j]);                              // pairs
      for (int i=0; i<agg.size(); i++)   
        for (int j=0; j<unagg.size(); j++)
          if (agg[i]==unagg[j]) {
            unagg.erase(unagg.begin()+j);                   // Remove any found from
            j--;                                            // unagg
          } 
    } 
  //analysis
    dist[agg.size()-1]++;                                     // One aggregate of size
    if (agg.size()+unagg.size()!=g.size()) {
      cout <<" Error in aggregation::count()"<<endl
           <<" g.size() = agg.size()+unagg.size()"<<endl
           <<"    "<< g.size()<<" = "<<agg.size()<<" +  "<<unagg.size()<<endl;
    }
    point CMagg, fix;
    for (int i=1; i<agg.size(); i++) {
      fix=agg[i]->cm-agg[0]->cm;
      CMagg+=fix;
    }
    CMagg+=agg[0]->cm;
    con->boundary(CMagg);
    for (int i=0; i<agg.size();i++) {
      for (int j=agg[i]->beg; j<agg[i]->end+1; j++) {
        RG2[agg.size()-1]+=con->sqdist(con->p[j], CMagg);
      }
    }
  }
}
void aggregation::write(string file) {
  ofstream f(file.c_str());
  if (f) {
    for (int i=0;i<g.size();i++) {
      if (dist[i]!=0)
        f << i+1 << "  "<<double(dist[i])/CNT/g.size() <<"  " <<RG2[i].avg() <<endl;
   
    }
    f.close();
  }
}

}//namespace
//
