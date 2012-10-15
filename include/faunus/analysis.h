#ifndef FAUNUS_ANALYSIS_H
#define FAUNUS_ANALYSIS_H

#ifndef SWIG
#include <faunus/common.h>
#include <faunus/average.h>
#include <faunus/physconst.h>
#include <faunus/group.h>
#include <faunus/space.h>
#include <faunus/point.h>
#endif

using std::vector;

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
     * This is the base class for analysis routines. Derived class must implement:
     * \li a descriptive name
     * \li _info()
     *
     * It is strongly recommended that derived classes implement:
     * \li a sample(...) function that uses run() to check if the analysis should be run or not.
     * \li the cite string to provide external information
     */
    class AnalysisBase {
      private:
        virtual string _info()=0; //!< info all classes must provide
        virtual void _test(UnitTest&);
      protected:
        char w;               //!< width of info
        unsigned long int cnt;//!< number of samples - increased for every run()==true.
        string name;          //!< descriptive name
        string cite;          //!< reference, url, doi etc. describing the analysis
        bool run();           //!< true if we should run, false of not (based on runfraction)
      public:
        AnalysisBase();
        virtual ~AnalysisBase();
        string info();       //!< Print info and results
        double runfraction;  //!< Chance that analysis should be run (default 1.0 = 100%)
        void test(UnitTest&);//!< Perform unit test
    };

    /*!
     * \brief General class for handling 2D tables - xy date, for example.
     * \author Mikael Lund
     * \date Lund 2011
     * \note Tx is used as the std::map key and which may be problematic due to direct floating
     *       point comparison (== operator). We have not experienced any issues with this, though.
     * \todo We get correct behavior, but is it really OK to have virtual functions in class templates??
     */
    template<typename Tx, typename Ty>
      class Table2D {
        protected:
          typedef std::map<Tx,Ty> Tmap;
          Ty count() {
            Ty cnt=0;
            for (auto &m : map)
              cnt+=m.second;
            return cnt;
          }
          Tx dx;
          Tmap map;
          string name;
        private:
          Tx round(Tx x) { return (x>=0) ? int( x/dx+0.5 )*dx : int( x/dx-0.5 )*dx; }
          virtual double get(Tx x) { return operator()(x); }
        public:
          enum type {HISTOGRAM, XYDATA};
          type tabletype;
          /*!
           * \brief Constructor
           * \param resolution Resolution of the x axis
           * \param key Table type: HISTOGRAM or XYDATA
           */
          Table2D(Tx resolution=0.2, type key=HISTOGRAM) {
            assert( resolution>0 );
            dx=resolution;
            name.clear();
            tabletype=key;
          }

          virtual ~Table2D() {}

          /*! \brief Access operator - returns reference to y(x) */
          Ty& operator() (Tx x) {
            return map[ round(x) ];
          }

          /*! \brief Save table to disk */
          void save(string filename) {
            if (tabletype==HISTOGRAM) {
              if (map.size()>0) map.begin()->second*=2;   // compensate for half bin width
              if (map.size()>1) (--map.end())->second*=2; // -//-
            }

            if (~map.empty()) {
              std::ofstream f(filename.c_str());
              f.precision(10);
              if (f) {
                f << "# Faunus 2D table: " << name << endl;
                for (auto m : map)
                  f << m.first << " " << get( m.first ) << endl;
              }
            }

            if (tabletype==HISTOGRAM) {
              if (map.size()>0) map.begin()->second/=2;   // restore half bin width
              if (map.size()>1) (--map.end())->second/=2; // -//-
            }
          }

          bool load(string filename) {
            std::ifstream f(filename.c_str());
            if (f) {
              map.clear();
              f.ignore(std::numeric_limits<std::streamsize>::max(),'\n'); // ignore first line
              while (!f.eof()) {
                Tx x;
                Ty y;
                f >> x >> y;
                operator()(x)=y;
              }
              return true;
            }
            return false;
          }
      };

    /*!
     * \brief Radial distribution analysis
     * \author Mikael Lund
     * \date Lund 2011
     *
     * This radial distribution is defined as \f$ g(r) = \rho(r) / \rho(\infty) \f$ where \f$ \rho \f$ are
     * the particle densities in spherical volume element \c rdr and in the bulk, respectively.
     *
     * Example:
     * \code
     * short cation = atom["Na"].id;
     * short anion = atom["Cl"].id;
     * Analysis::RadialDistribution<float,int> rdf(0.2); // 0.2 Å resolution
     * rdf.sample( myspace, mygroup, cation, anion );
     * rdf.save("rdf.dat");
     * \endcode
     */
    template<typename Tx=double, typename Ty=int>
      class RadialDistribution : public Table2D<Tx,Ty> {
        private:
          typedef Table2D<Tx,Ty> Ttable;
          virtual double volume(Tx x) {
            return 4./3.*pc::pi*( pow(x+0.5*this->dx,3) - pow(x-0.5*this->dx,3) );
          }
          double get(Tx x) {
            assert( volume(x)>0 );
            assert( this->count()>0 );
            if (bulkconc.cnt==0) bulkconc+=1;
            return (double)this->operator()(x) / volume(x) / (double)this->count() / bulkconc.avg()
              * this->map.size() * this->dx;
          }
          Average<double> bulkconc; //!< Average bulk concentration
        public:
          Tx maxdist; //!< Pairs with distances above this value will be skipped (default: infinity)

          /*!
           * \param res Resolution of X axis
           */
          RadialDistribution(Tx res=0.2) : Ttable(res,Ttable::HISTOGRAM) {
            this->name="Radial Distribution Function";
            maxdist=pc::infty;
          }
          /*!
           * \brief Sample radial distibution of two atom types
           * \param spc Simulation space
           * \param g Group to search
           * \param ida Atom id of first particle
           * \param idb Atom id of second particle
           */
          void sample(Space &spc, Group &g, short ida, short idb) {
            for (auto i=g.begin(); i!=g.end()-1; i++)
              for (auto j=i+1; j!=g.end(); j++)
                if ( (spc.p[*i].id==ida && spc.p[*j].id==idb) || (spc.p[*i].id==idb && spc.p[*j].id==ida) ) {
                  Tx r=spc.geo->dist(spc.p[*i], spc.p[*j]);
                  if (r<=maxdist)
                    this->operator() (r)++; 
                }
            double bulk=0;
            for (auto i : g)
              if (spc.p[i].id==ida || spc.p[i].id==idb)
                bulk++;
            bulkconc += bulk / spc.geo->getVolume();
          }
      };

    template<typename Tx=double, typename Ty=int>
      class LineDistribution : public RadialDistribution<Tx,Ty> {
        private:
          double volume(Tx x) { return 1; }
        public:
          LineDistribution(Tx res=0.2) : RadialDistribution<Tx,Ty>(res) {
            this->name="Line Distribution";
          }
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
        std::map< string, Average<double> > Rg2, Rg, Re2, Rs, Rs2, Rg2x, Rg2y, Rg2z;
        double gyrationRadiusSquared(const Group&, const Space &);
        Point vectorEnd2end(const Group&, const Space &);
        void _test(UnitTest&);
        string _info();
      public:
        PolymerShape();
        Point vectorgyrationRadiusSquared(const Group&, const Space &);
        void sample(const Group&, const Space&); //!< Sample properties of Group (identified by group name)
    };

    /*!
     * \brief Analyse charge multipoles and their fluctuations of groups
     * \author Anil Kurut
     * \date 2012
     *
     * This analysis class will analyse selected groups and calculate their net-charge, dipole moment as well
     * as their variances. It is possible to exclude certain atom types by added their names to an exclusionlist.
     * Several groups may be analysed - the sample() function will automatically identify different groups via
     * their names. The dipole moment is calculated with respect to the mass center.
     */
    class ChargeMultipole : public AnalysisBase {
      private:
        std::map< string, Average<double> > Z, Z2, mu, mu2;
        double charge(const Group&, const Space&);
        double dipole(const Group&, const Space&);
        virtual bool exclude(const particle&);  //!< Determines particle should be excluded from analysis
        string _info();
      public:
        ChargeMultipole();
        void sample(const Group&, const Space&); //!< Sample properties of Group (identified by group name)
        void sample(const vector<GroupMolecular>&, const Space&); //!< Sample properties of Group (identified by group name)
        std::set<string> exclusionlist; //!< Atom names listed here will be excluded from the analysis.
    };

    /*
       class VectorAlignment : public AnalysisBase {
       private:
       virtual Point convert(const Group&, const Space&); // Returns points calculated from group properties
       public:
       void sample(const Group&, const Group&, const Space&);
       };
       */

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
