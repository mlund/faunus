#ifndef FAUNUS_WIDOM_H
#define FAUNUS_WIDOM_H

#include <faunus/common.h>
#include <faunus/average.h>

namespace Faunus {
  class checkValue;
  class Space;

  namespace Analysis {

    class PolymerShape {
      private:
        std::map< string, Average<double> > Rg2, Rg, Re2;
        double gyrationRadiusSquared(const Group&, const Space &);
      public:
        void sample(const Group&, const Space &);
        string info();
    };

    /*! \brief Widom method for excess chemical potentials
     *  \author Mikael Lund
     *
     *  This class will use the ghost particle insertion technique
     *  to insert a collection of particles which, when summed, should
     *  have no net charge. This is used to calculate the mean excess
     *  chemical potential and activity coefficient.
     */
    class Widom {
      private:
        Space* spcPtr;
        Energy::Energybase* potPtr;
        Average<double> expsum; //!< Average of the excess chemical potential 
      protected:
        long unsigned int cnt;  //!< count test insertions
        p_vec g;                //!< List of ghost particles to insert (simultaneously)
      public:
        Widom(Space&, Energy::Energybase&);
        string info();                           //!< Print results of analysis
        void addGhost(particle);                 //!< Add particle to insert
        void addGhost(Space&);                   //!< All all species found in the container
        void sample(int=10);                     //!< Insert and analyse
        void check(UnitTest&);                   //!< Output checking
        double gamma();                          //!< Mean activity coefficient
        double muex();                           //!< Mean excess chemical potential
    };

    /*!
     * Single particle Widom insertion analysis including
     * charge re-scaling for electrostatics according to
     * Svensson and Woodward, Mol. Phys. 1988, 64(2), 247-259.
     * Currently, the inserted particle is a charged, hard sphere.
     *
     * \brief Single particle hard sphere Widom insertion with charge scaling
     * \author Martin Trulsson and Mikael Lund
     * \date Lund / Prague 2007-2008.
     * \note This is a direct conversion of the Widom routine found in the bulk.f
     *       program by Bolhuis/Jonsson/Akesson
     */
    class WidomScaled {
      private:
        p_vec g;         //!< list of test particles
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
        bool overlap(particle&, particle&, Space&); //!< Particle overlap test

      public:
        WidomScaled(int=10);        //!< Constructor, number of test insertions
        void add(particle);     //!< Add test particle
        void add(Space&);
        string info();          //!< Get results
        void insert(Space&, double=7.1); //!< Ghost insertion
        //void insert(container &, energybase &, vector<point> &); //!< Ghost insertion in a set of points
    };
  }//namespace
}//namespace
#endif
