#ifndef FAU_WIDOMMOD_H
#define FAU_WIDOMMOD_H

#include "faunus/energy.h"
#include "faunus/container.h"
#include "faunus/analysis.h"

namespace Faunus {
  /*! \brief Widom method for excess chemical potentials
   *  \author Mikael Lund
   *
   *  This class will use the ghost particle insertion technique
   *  to insert a collection of particles which, when summed, should
   *  have no net charge. This is used to calculate the mean excess
   *  chemical potential and activity coefficient.
   */
  class widom : public analysis {
    private:
      average<double> expsum; //!< Average of the excess chemical potential 
    protected:
      int ghostin;
      long long int cnt;          //< count test insertions
      vector<particle> g;  //!< List of ghost particles to insert (simultaneously)
    public:
      widom(int n=10) {
        ghostin=n;
        runfraction=1.0;
      }
      string info();                          //!< Print results of analysis
      void add(particle);                     //!< Add particle to insert
      void add(container &);                  //!< All all species found in the container
      void insert(container &, energybase &); //!< Insert and analyse
      double gamma() { return exp(muex()); }  //!< Mean activity coefficient
      double muex() { return -log(expsum.avg())/g.size(); } //!< Mean excess chemical potential
  };

  /*|
   * Single particle Widom insertion analysis including
   * charge re-scaling for electrostatics according to
   * Svensson and Woodward, Mol. Phys. 1988, 64(2), 247-259.
   * Currently, the inserted particle is a charged, hard sphere.
   *
   * \author Martin Trulsson and Mikael Lund
   * \date Lund / Prague 2007-2008.
   * \note This is a direct conversion of the Widom routine found in the bulk.f
   *       program by Bolhuis/Jonsson/Akesson
   */

  class widomSW : public analysis {
    private:
      vector<particle> g;         //!< list of test particles
      vector<double> chel;        //!< electrostatic
      vector<double> chhc;        //!< hard collision
      vector<double> chex;        //!< excess
      vector<double> chexw;       //!< excess
      vector<double> chtot;       //!< total
      vector< vector<double> > ewden;     //!< charging denominator
      vector< vector<double> > ewnom;     //!< charging nominator
      vector< vector<double> > chint;     //!< charging integrand
      vector<double> chid;                //!< ideal term
      vector<double> expuw;
      vector<int> ihc,irej;
      long long int cnt;          //< count test insertions
      int ghostin;                //< ghost insertions
      void init();
      bool overlap(particle &, particle &, container &); //!< Particle overlap test

    public:
      widomSW(int=10);        //!< Constructor, number of test insertions
      void add(particle);     //!< Add test particle
      void add(container &);
      string info();          //!< Get results
      void insert(container &, energybase &); //!< Ghost insertion
  };
}
#endif
