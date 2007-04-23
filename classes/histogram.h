#ifndef _HISTOGRAM_H
#define _HISTOGRAM_H

#include <cmath>
#include <string>
#include "point.h"
#include "average.h"
#include "xytable.h"

/*!
 * \brief Histogram class
 * \author Mikael Lund
 */
class histogram : private xytable<float,unsigned long int> {
  friend class rdf;
  private:
    unsigned long long int cnt;
  public:
    string comment; //!< User defined comment
    void show();    //!< Show results for all x

    //! \param res x value resolution
    //! \param min minimum x value
    //! \param max maximum x value
    histogram(float res=0.5, float min=0, float max=0)
      : xytable<float,unsigned long int>(res,min,max)
    { cnt=0; }
    
    /*! Example of histogram class
     * \example histogram-test.C
     */

    //! Increment bin for x value
    void add(float x) {
      (*this)(x)++;
      cnt++;
    }

    //! \brief Get bin for x value
    //! \return \f$ \frac{N(r)}{N_{tot}}\f$
    virtual float get(float x) { return (*this)(x)/float(cnt); }    
};
void histogram::show() {
  float g;
  for (float x=0; x<xmax(); x+=xres) {
    g=get(x);
    if (g!=0)
      cout << x << " " << g << "\n";
  };
}

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
    rdf(short, short, float=.5, float=0); 
    void update(vector<particle> &);      //!< Update histogram vector
    float get(float);                     //!< Get g(x)
};

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

/*!
 * Calculate all distances between between species 1
 * and 2 and update the histogram.
 */
void rdf::update(vector<particle> &p) {
  int k,n=p.size();
  for (int i=0; i<n-1; i++)
    for (int j=i+1; j<n; j++) 
      if ( (p[i].id==a && p[j].id==b)
          || (p[j].id==a && p[i].id==b) ) 
        (*this)( abs(p[i].dist(p[j]))/xres+0.5 )++;
}

float rdf::volume(float x) { return 4./3.*acos(-1)*( pow(x+xres,3)-pow(x,3) ); }

/*!
 *  Get g(x) from histogram according to
 *    \f$ g(x) = \frac{N(r)}{N_{tot}} \frac{ 3 } { 4\pi\left [ (x+xres)^3 - x^3 \right ] }\f$
 */
float rdf::get(float x) { return (*this)(x)/(cnt*volume(x)); }
#endif
