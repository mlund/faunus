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

namespace Faunus {
  class checkValue;
  class Space;

  /*!
   * \brief Namespace for analysis routines
   */
  namespace Analysis {

    /**
     * @brief Base class for analysis routines.
     *
     * This is the base class for analysis routines.
     * Derived class must implement:
     *
     * - a descriptive name
     * - `_info()`
     *
     * It is strongly recommended that derived classes implement:
     *
     * - a sample(...) function that uses `run()` to check if the
     *   analysis should be run or not.
     * - the `cite` string to provide external information
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

    /**
     * @brief General class for handling 2D tables - xy date, for example.
     * @date Lund 2011
     * @note `Tx` is used as the `std::map` key and which may be
     * problematic due to direct floating point comparison (== operator).
     * We have not experienced any issues with this, though.
     *
     * @todo We get correct behavior, but is it really OK to have
     * virtual functions in class templates??
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
          /**
           * @brief Constructor
           * @param resolution Resolution of the x axis
           * @param key Table type: HISTOGRAM or XYDATA
           */
          Table2D(Tx resolution=0.2, type key=XYDATA) {
            tabletype=key;
            setResolution(resolution);
          }

          void clear() { map.clear(); }

          void setResolution(Tx resolution) {
            assert( resolution>0 );
            dx=resolution;
            map.clear();
          }

          virtual ~Table2D() {}

          /** @brief Access operator - returns reference to y(x) */
          Ty& operator() (Tx x) {
            return map[ round(x) ];
          }

          /** @brief Save table to disk */
          void save(string filename) {
            if (tabletype==HISTOGRAM) {
              if (!map.empty()) map.begin()->second*=2;   // compensate for half bin width
              if (map.size()>1) (--map.end())->second*=2; // -//-
            }

            if (!map.empty()) {
              std::ofstream f(filename.c_str());
              f.precision(10);
              if (f) {
                f << "# Faunus 2D table: " << name << endl;
                for (auto m : map)
                  f << m.first << " " << get( m.first ) << endl;
              }
            }

            if (tabletype==HISTOGRAM) {
              if (!map.empty()) map.begin()->second/=2;   // restore half bin width
              if (map.size()>1) (--map.end())->second/=2; // -//-
            }
          }

          /*! Returns x at minumum y */
          Tx miny() {
            assert(!map.empty());
            Ty min=std::numeric_limits<Ty>::max();
            Tx x=0;
            for (auto &m : map)
              if (m.second<min) {
                min=m.second;
                x=m.first;
              }
            return x;
          }

          /*! Returns x at minumum y */
          Tx maxy() {
            assert(!map.empty());
            Ty max=std::numeric_limits<Ty>::min();
            Tx x=0;
            for (auto &m : map)
              if (m.second>max) {
                max=m.second;
                x=m.first;
              }
            return x;
          }
        
          /*! Returns x at minumum x */
          Tx minx() {
            assert(!map.empty());
            Tx x=0;
            for (auto &m : map) {
              x=m.first;
              break;
            }
            return x;
          }

          /**
           * @brief Load table from disk
           * @note The first line - used for comments - is ignored.
           * @todo Implement end bin compensation as in the save()
           * function when loading HISTOGRAMs
           */
          bool load(string filename) {
            std::ifstream f(filename.c_str());
            if (f) {
              map.clear();
              f.ignore(std::numeric_limits<std::streamsize>::max(),'\n'); // ignore first line
              while (!f.eof()) {
                Tx x;
                double y;
                f >> x >> y;
                operator()(x)=y;
              }
              return true;
            }
            return false;
          }
      };

    /*!
      \brief General class for penalty functions along a coordinate
      \date Malmo, 2011

      This class stores a penalty function, f(x), along a given coordinate, x,
      of type `Tcoordinate` which could be a distance, angle, volume etc.
      Initially f(x) is zero for all x.
      Each time the system visits x the update(x) function should be called
      so as to add the penalty energy, du. In the energy evaluation, the
      coordinate x should be associated with the extra energy f(x).
      This will eventually ensure uniform sampling.

      Example:

      ~~~
      PenaltyFunction<double> f(0.1,1000,6.0); // 0.1 kT penalty
      Point masscenter;           // some 3D coordinate...
      ...
      f.update(masscenter.z);     // update penalty energy for z component
      double u = f(masscenter.z); // get accumulated penalty at coordinate (kT)
      f.save("penalty.dat");      // save to disk
      ~~~

      In the above example, the penalty energy will be scaled by 0.5 if the
      sampling along the coordinate is less than 6 kT between the least and
      most likely position.
      This threshold check is carried out every 1000th call to `update()`.
      Note also that when the penalty energy is scaled, so is the threshold
      (also by a factor of 0.5).

      */
    template<typename Tcoord=float>
      class PenaltyFunction : public Table2D<Tcoord,double> {
        private:
          unsigned long long _cnt;
          int _Ncheck;
          double _kTthreshold;
          typedef Table2D<Tcoord,double> Tbase;
          typedef Table2D<Tcoord,unsigned long long int> Thist;
          Thist hist;
          Tcoord _du; //!< penalty energy
          std::string _log;
        public:
          /*!
           * \brief Constructor
           * \param penalty Penalty energy for each update (kT)
           * \param Ncheck Check histogram every Nscale'th step (put large number for no scaling, default)
           * \param kTthreshold Half penalty energy once this threshold in distribution has been reached
           * \param res Resolution of the penalty function (default 0.1)
           */
          PenaltyFunction(double penalty, int Ncheck=1e9, double kTthreshold=5, Tcoord res=0.1)
            : Tbase(res, Tbase::XYDATA), hist(res, Thist::HISTOGRAM) {
              Tbase::name="penalty";
              _cnt=0;
              _Ncheck=Ncheck;
              _kTthreshold=kTthreshold;
              _du=penalty;
              assert(Ncheck>0);
              _log="#   initial penalty energy = "+std::to_string(_du)+"\n";
            }

          /*! \brief Update penalty for coordinate */
          double update(Tcoord coordinate) {
            _cnt++;
            Tbase::operator()(coordinate)+=_du;  // penalize coordinate
            hist(coordinate)++;                  // increment internal histogram
            if ((_cnt%_Ncheck)==0) {             // if Ncheck'th time
              double deltakT=log( hist(hist.maxy()) / double(hist(hist.miny())) );
              assert(deltakT>0);
              std::ostringstream o;
              o << "#   n=" << _cnt << " dkT=" << deltakT;
              if (deltakT<_kTthreshold) {   // if histogram diff. is smaller than threshold
                _kTthreshold*=0.5;          // ...downscale threshold
                scale(0.5);                 // ...and penalty energy
                o << " update: du=" << _du << " threshold=" << _kTthreshold;
              }
              _log += o.str() + "\n";       // save info to log
            }
            return _du;
          }
          /*! \brief Manually scale penalty energy */
          void scale(double s) { _du*=s; }

          /*! \brief Save table to disk */
          void save(string filename) {
            Tbase::save(filename);
            hist.save(filename+".dist");
          }
        
          string info() {
            return "# Penalty function log:\n" + _log;
          }
      };

    template<typename Tx, typename Ty=unsigned long int>
      class Histogram : public Table2D<Tx,Ty> {
        public:
          Histogram(Tx resolution=0.2) : Table2D<Tx,Ty>(resolution, Table2D<Tx,Ty>::HISTOGRAM) {
            static_assert( std::is_integral<Ty>::value, "Histogram must be of integral type");
            static_assert( std::is_unsigned<Ty>::value, "Histogram must be unsigned");
          }
      };

    /*!
     * \brief Radial distribution analysis
     *
     * This radial distribution is defined as \f$ g(r) = \rho(r) / \rho(\infty) \f$ where \f$ \rho \f$ are
     * the particle densities in spherical volume element `rdr` and in the bulk, respectively.
     *
     * Example:
     *
     * ~~~
     * short cation = atom["Na"].id;
     * short anion = atom["Cl"].id;
     * Analysis::RadialDistribution<float,unsigned int> rdf(0.2); // 0.2 Å resolution
     * rdf.sample( myspace, mygroup, cation, anion );
     * rdf.save("rdf.dat");
     * ~~~
     *
     * \date Lund 2011
     */
    template<typename Tx=float, typename Ty=unsigned long long int>
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
            static_assert( std::is_integral<Ty>::value, "Histogram must be of integral type");
            static_assert( std::is_unsigned<Ty>::value, "Histogram must be unsigned");
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

    template<typename Tx=double, typename Ty=unsigned long>
      class LineDistribution : public RadialDistribution<Tx,Ty> {
        private:
          double volume(Tx x) { return 1; }
        public:
          LineDistribution(Tx res=0.2) : RadialDistribution<Tx,Ty>(res) {
            this->name="Line Distribution";
          }

      };

    /*!
     * \brief Line distr. when the bins values should sum up to `n`.
     *
     * Example: Salt line distribution
     *
     * ~~~
     * Analysis::LineDistributionNorm<float,unsigned long int> saltdistr(salt.size(), 0.2);
     * ~~~
     *
     * \author Axel Thuresson
     * \date Lund 2012
     */
    template<typename Tx=double, typename Ty=int>
    class LineDistributionNorm : public RadialDistribution<Tx,Ty> {
    private:
      double volume(Tx x) { return 1; }
      int n;
    public:
      LineDistributionNorm(int al_n=1, Tx res=0.2) : RadialDistribution<Tx,Ty>(res) {
        this->name="Line Distribution Normalized for n particles";
        n = al_n;
      }
      double get(Tx x) {
        assert( volume(x)>0 );
        assert( this->count()>0 );
        return (double)this->operator()(x) * n / (double)this->count();
      }
      
      /*!
       * \brief Simplest form of the midplane pressure
       */
      double mid() {
        return (get(this->dx)+get(-this->dx))*0.5/this->dx;
      }
      
      /*!
       * \brief Simplest form of the end pressure
       */
      double end() {
        return (get(this->minx())+get(this->minx()+this->dx)+get(-this->minx())+get(-this->minx()-this->dx))*0.25/this->dx;
      }
      
    };
    
    
    /*!
     * \brief Base class for force calculations
     *
     * Includes some neccessary functionality for deriving the force.
     *
     * \author Axel Thuresson
     * \date Lund, 2013
     */
    
    class TwobodyForce : public AnalysisBase {
    protected:
      Energy::Energybase* pot;         //!< Pointer to energy functions
      Space* spc;                      //!< Pointer to Space (particles and groups are stored there)
      string _info();         //!< Print results of analysis
      Group* igroup1;
      Group* igroup2;
      Group* ions;
    public:
      TwobodyForce(InputMap&, Energy::Energybase&, Space&, Group &, Group &, Group &);//!< Constructor
      virtual void calc();
      void save(string);
      void setTwobodies(Group &, Group &, Group &);
      virtual Point meanforce();
    };
    
    /*!
     * \brief Calculates the "direct" twobody mean force
     *
     * Calculates the "direct" mean force between two bodies including ions.
     * This method is usually decent at close distances (large mean force).
     * When the two bodies are far apart (small mean force) the difference 
     * between two large is taken which gives a relative large error.
     *
     * \author Axel Thuresson
     * \date Lund, 2013
     */
    
    class TwobodyForceDirect : public TwobodyForce {
    private:
      Point f_pp;
      Point f_pi;
      Point f_ip;
    protected:
      string _info();         //!< Print results of analysis
    public:
      TwobodyForceDirect(InputMap&, Energy::Energybase&, Space&, Group &, Group &, Group &);//!< Constructor
      void calc();
      Point meanforce();
    };
    
    
    
    /*!
     * \brief Calculates the midplane twobody mean force
     *
     * Calculates the midplane mean force between two bodies including ions.
     * This method has usually faster convergence than direct force calculations.
     *
     * \author Axel Thuresson
     * \date Lund, 2013
     */
    
    class TwobodyForceMidp : public TwobodyForce {
    private:
      Point f_pp;
      Point f_pi;
      Point f_ip;
      Point f_ii;
      Analysis::LineDistributionNorm<float,unsigned long int> *saltdistr;
    protected:
      string _info();         //!< Print results of analysis
    public:
      TwobodyForceMidp(InputMap&, Energy::Energybase&, Space&, Group &, Group &, Group &, Analysis::LineDistributionNorm<float,unsigned long int>*);//!< Constructor
      void calc();
      Point meanforce();
    };

    /*!
     * \brief Analysis of polymer shape - radius of gyration, shape factor etc.
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

    /**
     * @brief Analyse charge multipoles and their fluctuations of groups
     *
     * This analysis class will analyse selected groups and calculate
     * their net-charge, dipole moment as well as their variances.
     * It is possible to exclude certain atom types by added their
     * names to an exclusionlist. Several groups may be analysed -
     * the `sample()` function will automatically identify different
     * groups via their names.
     * The dipole moment is calculated with respect to the mass center.
     *
     * @author Anil Kurut
     * @date Lund 2012
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

        /* @brief Sample properties of Group (identified by group name) */
        template<typename Tgroup>
          void sample(const std::vector<Tgroup> &gvec, const Space &spc) {
            if (!run())
              return;
            for (auto &g : gvec)
              sample(g, spc);
          }

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

    /**
     * @brief Widom method for excess chemical potentials
     *
     * This class will use the ghost particle insertion technique
     * to insert a collection of particles which, when summed, should
     * have no net charge. This is used to calculate the mean excess
     * chemical potential and activity coefficient.
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
        void addGhost(particle);                 //!< Add particle to insert - sum of all added particle charges should be zero.
        void addGhost(Space&);                   //!< All species found in the container
        void sample(int=10);                     //!< Insert and analyse `n` times.
        void check(UnitTest&);                   //!< Output checking
        double gamma();                          //!< Sampled mean activity coefficient
        double muex();                           //!< Sampled mean excess chemical potential
    };

    /**
     * @brief Single particle hard sphere Widom insertion with charge scaling
     *
     * This will calculate excess chemical potentials for single particles
     * in the primitive model of electrolytes. Use the `add()` functions
     * to add test or *ghost* particles and call `insert()` to perform single
     * particle insertions.
     * Inserted particles are *non-perturbing* and thus removed again after
     * sampling. Electroneutrality for insertions of charged species is
     * maintaing by charge re-scaling according to 
     *
     * - [Svensson and Woodward, Mol. Phys. 1988, 64:247]
     *   (http://doi.org/ft9bv9)
     *
     * Currently this works **only** for the primitive model of electrolytes, i.e.
     * hard, charged spheres interacting with a Coulomb potential.
     *
     * @note This is a conversion of the Widom routine found in the `bulk.f`
     *       fortran program by Bolhuis/Jonsson/Akesson at Lund University.
     * @author Martin Trulsson and Mikael Lund
     * @date Lund / Prague 2007-2008.
     */
    class WidomScaled : public AnalysisBase {
      private:
        typedef std::vector<double> Tvec;
        string _info();     //!< Get results
        p_vec g;            //!< list of test particles
        Tvec chel;          //!< electrostatic
        Tvec chhc;          //!< hard collision
        Tvec chex;          //!< excess
        Tvec chexw;         //!< excess
        Tvec chtot;         //!< total
        vector<Tvec> ewden; //!< charging denominator
        vector<Tvec> ewnom; //!< charging nominator
        vector<Tvec> chint; //!< charging integrand
        Tvec chid;          //!< ideal term
        Tvec expuw;
        vector<int> ihc,irej;
        long long int cnt;  //< count test insertions
        int ghostin;        //< ghost insertions
        void init();
        bool overlap(const particle&, const particle&, const Geometry::Geometrybase&); //!< Test overlap
        double lB;          //!< Bjerrum length [a]
      public:
        WidomScaled(double,int=10); //!< Constructor
        void add(const particle&);  //!< Add test particle type
        void add(const p_vec&);     //!< Add all unique particle types present in vector

        /** @brief Do test insertions and sample excess chemical potential */
        void insert(const p_vec&, Geometry::Geometrybase&);
    };
  }//namespace
}//namespace
#endif
