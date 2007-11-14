#ifndef _HISTOGRAM_H
#define _HISTOGRAM_H

#include <cmath>
#include <string>
#include <iostream>
#include "point.h"
#include "average.h"
#include "xytable.h"
#include "container.h"

//---------------------------------------------------
/*!
 * \brief Histogram class
 * \author Mikael Lund
 */
class histogram : private xytable<float,unsigned long int> {
  friend class rdf;
  friend class rdfP3;
  private:
    unsigned long int cnt;
  public:
    histogram(float, float, float);
    string comment;                     //!< User defined comment
    void add(float); 
    void write(string); 
    virtual float get(float); 

    //! Example of histogram class
    //! \example histogram-test.C
};
//! \param res x value resolution
//! \param min minimum x value
//! \param max maximum x value
histogram::histogram(float res, float min, float max)
: xytable<float,unsigned long int>(res,min,max) { cnt=0; }

//! Increment bin for x value
void histogram::add(float x) {
  (*this)(x)++;
  cnt++;
}

//! Get bin for x value
//! \return \f$ \frac{N(r)}{N_{tot}}\f$
float histogram::get(float x) { return (*this)(x)/float(cnt); }

//! Show results for all x
void histogram::write(string file) {
  float g;
  ofstream f(file.c_str());
  if (f) {
    for (float x=0; x<xmax(); x+=xres) {
      g=get(x);
      if (g!=0)
        f << x << " " << g << "\n";
    }
    f.close();
  }
}

//-----------------------------------------------------------------
/*!
 *  \brief Class to calculate the radial distribution between particles.
 *
 *  The internal histogram vector will be automatically expanded but an initial
 *  maximum x valued can be specified so as to utilize memory more efficiently.
 *  \author Mikael Lund
 *  \date Prague, April 2007
 *  \todo Needs testing!
 *
 */
class rdf : public histogram {
  private:
    short a,b;                   //!< Particle types to investigate
    float volume(float);         //!< Volume of shell r->r+xres
  public:
    rdf(float=.5, float=0);
    rdf(short, short, float=.5, float=0); 
    void update(container &);             //!< Update histogram vector
    void update(vector<particle> &);      //!< Update histogram vector
    void update(container &, point &, point &); //!< Update for two points
    float get(float);                     //!< Get g(x)
};

/*
 * \brief RDF class for cubic min. image
 * \todo Obselete?
 */
class rdfP3 : public histogram {
  private:
    short a,b;                   //!< Particle types to investigate
    float volume(float);         //!< Volume of shell r->r+xres
  public:
    rdfP3(short, short, float=.5, float=0, float=0);
    void update(vector<particle> &);      //!< Update histogram vector
    void update(vector<particle*> &);
    float get(float);                     //!< Get g(x)
    double len, len_inv;
};

/*!
 * \param resolution Histogram resolution (binwidth)
 * \param xmaximum Maximum x value (used for better memory utilisation)
 */
rdf::rdf( float resolution, float xmaximum) :
  histogram(resolution, 0, xmaximum) { }

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
rdfP3::rdfP3(short species1, short species2, float resolution, float xmaximum, float boxl) :
  histogram(resolution, 0, xmaximum)
{
  a=species1;
  b=species2;
  len=boxl;
  len_inv=1/boxl;
}

/*!
 * Update histogram between two known points
 *
 * \note Uses the container distance function
 */
void rdf::update(container &c, point &a, point &b) {
  add( abs( c.dist(a, b) ) );
}

/*!
 * Calculate all distances between between species 1
 * and 2 and update the histogram.
 *
 * \note Uses the container function to calculate distances
 */
void rdf::update(container &c) {
  unsigned short i,j,n=c.p.size();
  for (i=0; i<n-1; i++)
    for (j=i+1; j<n; j++) 
      if ( (c.p[i].id==a && c.p[j].id==b)
          || (c.p[j].id==a && c.p[i].id==b) )
        update( c, c.p[i], c.p[j] );
}

/*!
 * Calculate all distances between between species 1
 * and 2 and update the histogram.
 *
 * \warning This function uses a simple distance function (no mim. image)
 */
void rdf::update(vector<particle> &p)
{
  unsigned short i,j,n=p.size();
  for (i=0; i<n-1; i++)
    for (j=i+1; j<n; j++) 
      if ( (p[i].id==a && p[j].id==b)
          || (p[j].id==a && p[i].id==b) )
        add( abs(p[i].dist(p[j])) );
}

/*!
 * Calculate all distances between vector<particle> and 
 * update the histogram under P3 conditions
 */
void rdfP3::update(vector<particle*> &p)
{
  unsigned short i,j,n=p.size();
  for (i=0; i<n-1; i++)
    for (j=i+1; j<n; j++) 
      add( abs(p[i]->dist(*p[j], len, len_inv )));
}
float rdf::volume(float x) { return 4./3.*acos(-1)*( pow(x+xres,3)-pow(x,3) ); }
/*!
 *  Get g(x) from histogram according to
 *    \f$ g(x) = \frac{N(r)}{N_{tot}} \frac{ 3 } { 4\pi\left [ (x+xres)^3 - x^3 \right ] }\f$
 */
float rdf::get(float x) { return (*this)(x)/(cnt*volume(x)); }

float rdfP3::volume(float x) { return 4./3.*acos(-1)*( pow(x+xres,3)-pow(x,3) ); }
/*!
 *  Get g(x) from histogram according to
 *    \f$ g(x) = \frac{N(r)}{N_{tot}} \frac{ 3 } { 4\pi\left [ (x+xres)^3 - x^3 \right ] }\f$
 */
float rdfP3::get(float x) { return (*this)(x)/(cnt*volume(x)); }
#endif
