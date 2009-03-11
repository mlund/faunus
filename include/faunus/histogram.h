#ifndef FAU_HISTOGRAM_H
#define FAU_HISTOGRAM_H

#include "faunus/common.h"
#include "faunus/point.h"
#include "faunus/average.h"
#include "faunus/xytable.h"
#include "faunus/container.h"
#include "faunus/group.h"

namespace Faunus {
  /*!
   * \brief Histogram class
   * \author Mikael Lund
   */
  class histogram : private xytable<float,unsigned long int> {
    friend class FAUrdf;
    private:
    unsigned long int cnt;
    float xmaxi,xmini;  // ugly solution!
    public:
    histogram(float, float, float);
    string comment;                     //!< User defined comment
    void add(float);
    void write(string);
    virtual float get(float);

    //! Example of histogram class
    //! \example histogram-test.C
  };

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
  class FAUrdf : public histogram {
    private:
      short a,b;                   //!< Particle types to investigate
      float volume(float);         //!< Volume of shell r->r+xres
      unsigned int npart;
    public:
      FAUrdf(float=.5, float=0);
      FAUrdf(short, short, float=.5, float=0); 
      void update(container &, point &, string);//!< Search around a point
      void update(container &);                 //!< Update histogram vector
      void update(container &, group &);        //!< Update histogram vector - search only in group
      void update(vector<particle> &);         //!< Update histogram vector
      void update(container &, point &, point &);//!< Update for two points
      float get(float);                        //!< Get g(x)
  };
};//namespace
#endif
