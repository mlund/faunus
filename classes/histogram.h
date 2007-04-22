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
    histogram(float res=0.5, float min=0, float max=0)
      : xytable<float,unsigned long int>(res,min,max)
    { cnt=0; }
    
    //! Increase count for bin 
    void add(float r) {
      (*this)(r)++;
      cnt++;
    }

    //! Get value of bin (divided by the total counts)
    float get(float r) { return (*this)(r)/float(cnt); }
};

/*!
 *  \brief Class to calculate the radial distribution function between particles.
 *  \author Mikael Lund
 *  \date Prague, April 2007
 *  \todo Needs testing. Could instead be based on the xydata class.
 */
class rdf : public histogram {
  private:
    short a,b;                  //!< Particle types to investigate
    float volume(float);         //!< Volume of shell r->r+xres
  public:
    rdf(short, short, float=.5, float=0);      //!< Constructor
    void show();                                //!< Print g(r) for all r
    void update(vector<particle> &);            //!< Update histogram vector
    float get(float);                        //!< Get g(r)
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

/*!
 *  Get g(r) from histogram according to
 *    \f$ g(r) = \frac{N(r)}{N_{tot} 4\pi/3 \left ( (R_{cell}+binwidth)^3 - R_{cell}^3 \right ) } \f$
 */
float rdf::volume(float r) { return 4./3.*acos(-1)*( pow(r+xres,3)-pow(r,3) ); }
float rdf::get(float r) { return (*this)(r)/(cnt*volume(r)); }
void rdf::show() {
  float g;
  for (float r=0; r<xmax(); r+=xres) {
    g=get(r);
    if (g!=0)
      cout << r << " " << g << "\n";
  };
}

#endif
