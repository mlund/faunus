#ifndef _RDF_H
#define _RDF_H

#include <vector>
#include <cmath>
#include "point.h"

/*!
 *  \brief Class to calculate the radial distribution function between particles.
 *  \author Mikael Lund
 *  \date Prague, April 2007
 *  \todo Needs testing. Could instead be based on the xydata class.
 */
class rdf {
  private:
    short a,b;                  //!< Particle types to investigate
    double res;                 //!< Width of histogram bins (resolution).
    double bulkc;               //!< Concentration in the bulk
    unsigned long long int cnt; //!< Total number of counts in the histrogram
    vector<unsigned long long int> v;           //!< Histogram vector
  public:
    rdf(short, short, double=.5, double=1);     //!< Constructor
    void show();                                //!< Print g(r) for all r
    void update(vector<particle> &);            //!< Update histogram vector
    double gofr(double);                        //!< Get g(r)
};

/*!
 * \param species1 Particle type 1
 * \param species2 Particle type 2
 * \param resolution Histogram resolution (binwidth)
 * \param bulkconc Concentration in the bulk
 */
rdf::rdf(short species1, short species2, double resolution, double bulkconc) {
  cnt=0;
  a=species1;
  b=species2;
  res=resolution;
  bulkc=bulkconc;
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
          || (p[j].id==a && p[i].id==b) )  {
        k=static_cast<int>(p[i].dist(p[j])/res+0.5);
        if (k>=v.size())
          v.resize(k+1);
        v[k]++;
        cnt++;
      };
}

/*!
 *  Get g(r) from histogram according to
 *    \f$ g(r) = \frac{N(r)}{N_{tot} 4\pi/3 \left ( (R_{cell}+binwidth)^3 - R_{cell}^3 \right ) } \f$
 */
double rdf::gofr(double r) {
  int k=static_cast<int>(r/res+0.5);
  double vol=4./3.*acos(-1)*( pow( (k+1)*res,3 )-pow( k*res,3 ) );
  return static_cast<double>(v[k]) / cnt / vol;
}

void rdf::show() {
  double g;
  for (int i=0; i<v.size(); i++) {
    g=gofr(i*res)/bulkc;
    if (g!=0)
      cout << i*res << " " << g << "\n";
  };
}
#endif
