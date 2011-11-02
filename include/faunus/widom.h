#ifndef FAUNUS_WIDOM_H
#define FAUNUS_WIDOM_H

#include <faunus/common.h>
#include <faunus/average.h>

namespace Faunus {
  class checkValue;
  class Space;

  /*!
   * \brief Namespace for analysis routines
   */
  namespace Analysis {

    /*!
     * \brief Base class for analysis routines.
     *
     * Derived classes are expected to provide a name and _info()
     * function. It is recommended that derived classes also implement
     * a sample(...) function that uses the run() function to check if the
     * analysis should be run or not.
     */
    class AnalysisBase {
      private:
        virtual string _info()=0; //!< info all classes must provide
      protected:
        char w;               //!< width of info
        unsigned long int cnt;//!< number of samples - increased for every run()==true.
        string name;          //!< descriptive name
        string cite;          //!< reference, url, doi etc. describing the analysis
        bool run();           //!< true if we should run, false of not (based on runfraction)
      public:
        AnalysisBase();
        string info();       //!< Print info and results
        double runfraction;  //!< Chance that analysis should be run (default 1.0 = 100%)
    };

    /*!
     * \brief Analysis of polymer shape - radius of gyration, shape factor etc.
     * \author Mikael Lund
     * \date November, 2011
     *
     * This will analyse polymer Groups and calculate Rg, Re and the shape factor. If
     * sample() is called with different groups these will be distinguished by their
     * *name* and sampled individually.
     */
    class PolymerShape : public AnalysisBase {
      private:
        std::map< string, Average<double> > Rg2, Rg, Re2;
        double gyrationRadiusSquared(const Group&, const Space &);
        string _info();
      public:
        PolymerShape();
        void sample(const Group&, const Space&); //!< Sample properties of Group (identified by group name)
    };

    /*! \brief Widom method for excess chemical potentials
     *  \author Mikael Lund
     *
     *  This class will use the ghost particle insertion technique
     *  to insert a collection of particles which, when summed, should
     *  have no net charge. This is used to calculate the mean excess
     *  chemical potential and activity coefficient.
     */
    class Widom : public AnalysisBase {
      private:
        Space* spcPtr;
        Energy::Energybase* potPtr;
        Average<double> expsum; //!< Average of the excess chemical potential 
        string _info();         //!< Print results of analysis
      protected:
        p_vec g;                //!< List of ghost particles to insert (simultaneously)
      public:
        Widom(Space&, Energy::Energybase&);
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
    class WidomScaled : public AnalysisBase {
      private:
        string _info();   //!< Get results
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
        void insert(Space&, double=7.1); //!< Ghost insertion
        //void insert(container &, energybase &, vector<point> &); //!< Ghost insertion in a set of points
    };
  }//namespace
}//namespace
#endif
