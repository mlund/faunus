#include <faunus/histogram.h>
#include <faunus/io.h>
#include <faunus/species.h>
#include <faunus/point.h>
#include <faunus/geometry.h>
#include <faunus/group.h>
#include <faunus/space.h>

namespace Faunus {
  //! \param res x value resolution
  //! \param min minimum x value
  //! \param max maximum x value
  histogram::histogram(float res, float min, float max)
    : xytable<float,unsigned long int>(res,min,max) {
      reset(res,min,max);
    }
  
  void histogram::reset(float res, float min, float max) {
    cnt=0;
    xmaxi=max;
    xmini=min;
    init(res,min,max);
  }

  //! Increment bin for x value
  void histogram::add(float x) {
    if (x>=xmaxi || x<=xmini) return;
    (*this)(x)++;
    cnt++;
  }

  //! Get bin for x value
  //! \return \f$ \frac{N(r)}{N_{tot}}\f$
  float histogram::get(float x) {
    return (*this)(x)/float(cnt);
  }

  vector<double> histogram::xvec() {
    vector<double> v;
    for (float x=xmin; x<xmax(); x+=xres)
      v.push_back(x);
    return v;
  }

  vector<double> histogram::yvec() {
    vector<double> v;
    for (float x=xmin; x<xmax(); x+=xres)
      v.push_back( get(x) );
    return v;
  }

  //! Show results for all x
  void histogram::write(string file) {
    float g;
    std::ofstream f(file.c_str());
    if (f) {
      f.precision(6);
      for (double x=xmin; x<xmax(); x+=xres) {
        g=get(x);
        if (g!=0.0) {
          if (x+xres>=xmax() || x==xmin) // double the very
            g=g*2;                       // first and last points
          f << x << " " << g << "\n";
        }
      }
      f.close();
    }
  }

  //! Dump table to disk
  void histogram::save(string filename) {
    xytable<float,unsigned long int>::save(filename);
  }

  /*!
   * \param species1 Particle type 1
   * \param species2 Particle type 2
   * \param resolution Histogram resolution (binwidth)
   * \param xmaximum Maximum x value (used for better memory utilisation)
   */
  rdf::rdf(short species1, short species2, float resolution, float xmaximum) :
    histogram(resolution, 0, xmaximum)
  {
    a=species1;
    b=species2;
  }
  rdf::rdf(float resolution, float xmaximum, float xminimum) : histogram(resolution, xminimum, xmaximum) {}

  /*!
   * Update histogram between two known points
   * \note Uses the space distance function
   */
  void rdf::update(space &s, point &a, point &b) {
    add( s.geo->_dist(a, b) );
  }

  /*!
   * Calculate all distances between between species 1
   * and 2 and update the histogram.
   *
   * \note Uses the container function to calculate distances
   */
  void rdf::update(space &c) {
    int n=c.p.size();
    npart=0;
#pragma omp parallel for schedule (dynamic)
    for (int i=0; i<n-1; i++)
      for (int j=i+1; j<n; j++) 
        if ( (c.p[i].id==a && c.p[j].id==b)
            || (c.p[j].id==a && c.p[i].id==b) ) {
          npart++;
          update( c, c.p[i], c.p[j] );
        }
  }

  void rdf::update(space &c, group &g) {
    int n=g.end+1;
    npart=0;
    //#pragma omp for
    for (int i=g.beg; i<n-1; i++) {
      for (int j=i+1; j<n; j++) { 
        if ( (c.p[i].id==a && c.p[j].id==b)
            || (c.p[j].id==a && c.p[i].id==b) ) {
          npart++;
          add( c.geo->_dist(c.p[i], c.p[j]) );
        }
      }
    }
  }

  void rdf::update(space &c, point &p, string name) {
    int id=atom[name].id, n=c.p.size();
    npart=0;
    //#pragma omp for
    for (int i=0; i<n; ++i)
      if ( c.p[i].id==id ) {
        npart++;
        add( c.geo->_dist(p, c.p[i] ));
      }
  }

  /*!
   * Calculate shell volume at x
   */
  float rdf::volume(float x) { return 4./3.*acos(-1.)*( pow(x+0.5*xres,3)-pow(x-0.5*xres,3) ); }
  /*!
   *  Get g(x) from histogram according to
   *    \f$ g(x) = \frac{N(r)}{N_{tot}} \frac{ 3 } { 4\pi\left [ (x+xres)^3 - x^3 \right ] }\f$
   */
  float rdf::get(float x) {
    return (*this)(x)/(cnt*volume(x)) * 1660.57;
  }

  /*!
   *  Get g(x) from histogram according to
   *    \f$ g(x) = \frac{N(r)}{N_{tot}} \frac{ 3 } { 4\pi\left [ (x+xres)^3 - x^3 \right ] }\f$
   */

  profile::profile(float min, float max, float res) :
    xytable<float,unsigned long int>(res,min,max) { cnt=0; }
  float profile::conc(float x) { return ((*this)(x)>0) ? (*this)(x)/(cnt*volume(x)) : 0; }
  void profile::update(p_vec &p) { for (int i=0; i<p.size(); i++) add(p[i]); }
  bool profile::write(string name) {
    io fio;
    std::ostringstream o;
    for (float d=xmin; d<xmax(); d+=xres)
      o << d << " " << conc(d) << endl;
    return fio.writefile(name,o.str());
  }

  float cummsum::volume(double z) { return 1.; }

  cummsum::cummsum(
      unsigned char type, particle &center, float max,float res) :
    profile(0,max,res) {
      id=type;
      origo=&center;
    }

  //void cummsum::add(particle &p) { if (p.id==id) (*this)( sqrt(origo->sqdist(p)) )++; }

  float cylindric_profile::volume(double z) { return xres*acos(-1.)*r*r; }

  cylindric_profile::cylindric_profile(
      float radius, unsigned char type, float min,float max,float res) :
    profile(min,max,res) {
      r=radius;
      id=type;
    }

  void cylindric_profile::add(particle &p) {
    if (p.id==id)
      if (p.z>=xmin && p.z<=xmax())
        if (p.x*p.x+p.y*p.y<=r*r) {
          (*this)(p.z)++;
          cnt++;
        }
  }
  radial_profile::radial_profile(float min, float max, float res) :
    xytable<float,unsigned long int>(res,min,max) { cnt=0; }
  float radial_profile::volume(float x) {return acos(-1.)*(pow(x+xres*0.5,2.)-pow(x-xres*0.5,2.));}
  float radial_profile::conc(float x) { return ((*this)(x)>0) ? (*this)(x)/(cnt*volume(x)) : 0; }
  void  radial_profile::add(particle &p) { if (p.id==id){cnt++, (*this)(sqrt(pow(p.x-origo.x,2)+pow(p.y-origo.y,2)))++;} }
  void  radial_profile::add(point &o, point &p) { 
    cnt++ ;
    (*this)(sqrt(pow(p.x-o.x,2)+pow(p.y-o.y,2)))++; 
  }
  void  radial_profile::update(p_vec &p) { for (int i=0; i<p.size(); i++) add(p[i]); }
  bool  radial_profile::write(string name) {
    io fio;
    std::ostringstream o;
    for (float d=xmin; d<xmax(); d+=xres)
      o << d << " " << conc(d) << endl;
    return fio.writefile(name,o.str());
  }

  atomicRdf::atomicRdf(float dx, float max) : histogram(dx, 0, max) { }

  /*
  void atomicRdf::update(p_vec &p, group &g1, group &g2) {
    for (int i=g1.beg; i<=g1.end; i++)
      for (int j=g2.beg; j<=g2.end; j++)
        add( sqrt(p[i].sqdist(p[j])) );
  }
  */
}//namespace
