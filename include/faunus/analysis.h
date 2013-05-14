#ifndef FAUNUS_ANALYSIS_H
#define FAUNUS_ANALYSIS_H

#ifndef SWIG
#include <faunus/common.h>
#include <faunus/average.h>
#include <faunus/physconst.h>
#include <faunus/group.h>
#include <faunus/space.h>
#include <faunus/point.h>
#include <faunus/textio.h>
#endif

namespace Faunus {
  class checkValue;
  class Space;

  /**
   * @brief Namespace for analysis routines
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
          bool load(const string &filename) {
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

    /**
      @brief General class for penalty functions along a coordinate
      @date Malmo, 2011

      This class stores a penalty function, f(x), along a given coordinate, x,
      of type `Tcoordinate` which could be a distance, angle, volume etc.
      Initially f(x) is zero for all x.
      Each time the system visits x the update(x) function should be called
      so as to add the penalty energy, du. In the energy evaluation, the
      coordinate x should be associated with the extra energy f(x).

      This will eventually ensure uniform sampling. Example:

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
          /**
           * @brief Constructor
           * @param penalty Penalty energy for each update (kT)
           * @param Ncheck Check histogram every Nscale'th step
           *        (put large number for no scaling, default)
           * @param kTthreshold Half penalty energy once this
           *        threshold in distribution has been reached
           * @param res Resolution of the penalty function (default 0.1)
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

          /** @brief Update penalty for coordinate */
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
          void save(const string &filename) {
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

          /**
           * @brief Sample radial distibution of two atom types
           * @param spc Simulation space
           * @param g Group to search
           * @param ida Atom id of first particle
           * @param idb Atom id of second particle
           */
          template<class Tspace, class Tgroup>
            void sample(Tspace &spc, Tgroup &g, short ida, short idb) {
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
        void _test(UnitTest&);
        string _info();
        template<class Tgroup, class Tspace>
          double gyrationRadiusSquared(const Tgroup &pol, const Tspace &spc) {
            assert( spc.geo->dist(pol.cm, pol.massCenter(spc))<1e-9
                && "Mass center must be in sync.");
            Point rg2=vectorgyrationRadiusSquared(pol,spc);
            return rg2.x()+rg2.y()+rg2.z();
          }

        template<class Tgroup, class Tspace>
          Point vectorEnd2end(const Tgroup &pol, const Tspace &spc) {
            return spc.geo->vdist( spc.p[pol.front()], spc.p[pol.back()] );
          }

        template<class Tgroup, class Tspace>
          Point vectorgyrationRadiusSquared(const Tgroup &pol, const Tspace &spc) {
            assert( spc.geo->dist(pol.cm, pol.massCenter(spc))<1e-9
                && "Mass center must be in sync.");
            double sum=0;
            Point t, r2(0,0,0);
            for (auto i : pol) {
              t = spc.p[i]-pol.cm;                // vector to center of mass
              spc.geo->boundary(t);               // periodic boundary (if any)
              r2.x() += spc.p[i].mw * t.x() * t.x();
              r2.y() += spc.p[i].mw * t.y() * t.y();
              r2.z() += spc.p[i].mw * t.z() * t.z();
              sum += spc.p[i].mw;                 // total mass
            }
            assert(sum>0 && "Zero molecular weight not allowed.");
            return r2*(1./sum);
          }

      public:
        PolymerShape();

        /** @brief Sample properties of Group (identified by group name) */
        template<class Tgroup, class Tspace>
          void sample(const Tgroup &pol, const Tspace &spc) {
            if (!run() || pol.front()==pol.back())
              return;
            Point r2 = vectorgyrationRadiusSquared(pol,spc);
            double rg2 = r2.x()+r2.y()+r2.z(); 
            double re2 = spc.geo->sqdist( spc.p[pol.front()], spc.p[pol.back()] );
            Rg2[pol.name]  += rg2;
            Rg2x[pol.name] += r2.x();
            Rg2y[pol.name] += r2.y();
            Rg2z[pol.name] += r2.z();
            Rg[pol.name]   += sqrt(r2.x()+r2.y()+r2.z());
            Re2[pol.name]  += re2; //end-2-end squared
            double rs = Re2[pol.name].avg()/Rg2[pol.name].avg(); // fluctuations in shape factor
            Rs[pol.name]   += rs;
            Rs2[pol.name]  += rs*rs;
            //Point re = vectorEnd2end(pol,spc);
            //Re2[pol.name] += pow(re.len(), 2);
          }
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

        template<class Tgroup, class Tpvec>
          double charge(const Tgroup &g, const Tpvec &p, double Z=0) {
            for (auto i : g)
              if (!excluded(p[i]))
                Z+=p[i].charge;
            return Z;
          }

        template<class Tgroup, class Tspace>
          double dipole(const Tgroup &g, const Tspace &spc) {
            assert( spc.geo->dist(g.cm, g.massCenter(spc))<1e-9
                && "Mass center must be in sync.");
            Point t, mu(0,0,0);
            for (auto i : g)
              if (!excluded(spc.p[i])) {
                t = spc.p[i]-g.cm;                // vector to center of mass
                spc.geo->boundary(t);               // periodic boundary (if any)
                mu+=spc.p[i].charge*t;
              }
            return mu.len();
          }

        /** @brief Determines particle should be excluded from analysis */
        template<class Tparticle>
          bool excluded(const Tparticle &p){
            if (exclusionlist.count(atom[p.id].name)==0)
              return false;
            return true;
          }

        string _info();
      public:
        ChargeMultipole();

        /** @brief Sample properties of Group (identified by group name) */
        template<class Tgroup, class Tspace>
          void sample(const Tgroup &g, const Tspace &spc) {
            assert(!g.name.empty() && "All Groups should have a name!");
            if (run()) {
              double z=charge(g, spc.p);
              Z[g.name]+=z;
              Z2[g.name]+=pow(z,2);
              double dip=dipole(g,spc);
              mu[g.name]+=dip;
              mu2[g.name]+=pow(dip,2);
            }
          }

        /* @brief Sample properties of Group (identified by group name) */
        template<typename Tgroup, typename Tspace>
          void sample(const std::vector<Tgroup> &gvec, const Tspace &spc) {
            if (run())
              for (auto &g : gvec)
                sample(g, spc);
          }

        std::set<string> exclusionlist; //!< Atom names listed here will be excluded from the analysis.
    };

    /**
     * @brief Widom method for excess chemical potentials
     *
     * This class will use the ghost particle insertion technique
     * to insert a collection of particles which, when summed, should
     * have no net charge. This is used to calculate the mean excess
     * chemical potential and activity coefficient.
     */
    template<class Tparticle>
      class Widom : public AnalysisBase {
        private:
          Average<double> expsum; //!< Average of the excess chemical potential 

          string _info() {
            using namespace Faunus::textio;
            std::ostringstream o;
            o << pad(SUB,w, "Number of insertions") << expsum.cnt << endl
              << pad(SUB,w, "Excess chemical pot.") << muex() << kT << endl
              << pad(SUB,w, "Mean activity coefficient") << gamma() << endl
              << pad(SUB,w, "Ghost particles");
            for (auto &p : g)
              o << atom[p.id].name << " ";
            return o.str() + "\n";
          }

          void _test(UnitTest &test) { test("widom_muex", muex() ); }

        protected:
          std::vector<Tparticle> g; //!< Pool of ghost particles to insert (simultaneously)
        public:
          Widom() {
            name="Multi Particle Widom Analysis";
            cite="doi:10/dkv4s6";
          }

          void addGhost(Tparticle p) { g.push_back(p); }

          /* @brief Add particle to insert - sum of added particle charges should be zero.*/
          template<class Tpvec>
            void addGhost(Tpvec &p) {
              std::map<short,bool> map; // replace w. `std::set`
              for (auto i : p)
                map[ i.id ] = true;
              for (auto &m : map) {
                particle a;
                a=atom[m.first];
                addGhost(a);
              }
            }

          /** @brief Sampled mean activity coefficient */
          double gamma() { return exp(muex()); }

          /** @brief Sampled mean excess chemical potential */
          double muex() { return -log(expsum.avg())/g.size(); }

          /** @brief Insert and analyse `n` times */
          template<class Tspace, class Tenergy>
            void sample(int ghostin, Tspace &spc, Tenergy &pot) {
              if (!run())
                return;
              assert(spc.geo!=NULL);
              int n=g.size();
              for (int k=0; k<ghostin; k++) {     // insert ghostin times
                double du=0;
                for (int i=0; i<n; i++)
                  spc.geo->randompos( g[i] ); // random ghost positions
                for (int i=0; i<n; i++)
                  pot.all2p( spc.p, g[i] );    // energy with all particles in space
                for (int i=0; i<n-1; i++)
                  for (int j=i+1; j<n; j++)
                    du+=pot.p2p( g[i], g[j] );   // energy between ghost particles
                expsum += exp(-du);
              }
            }
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

    /**
     * @brief Samples bilayer structure
     *
     * This was developed for coarse grained membrane models
     * but should be general enough for other uses.
     */
    class BilayerStructure : public AnalysisBase {

      private:

        inline string _info() {
          using namespace textio;
          std::ostringstream o;
          if (cnt>0)
            o << pad(SUB,w,"Lipid order parameter") << S << endl
              << pad(SUB,w,"Area per lipid") << A << " "+sigma+squared << endl;
          return o.str();
        }

        Average<double> S, A;

        void _test(UnitTest &t);

      public:

        inline BilayerStructure() {
          name="Bilayer structure";
          cite="doi:10/chqzjk";
        }

        template<class Tcuboid, class Tpvec, class Tgroup>
          void sample(Tcuboid &geo, Tpvec &p, Tgroup &lipids) {
            if (run()) {
              cnt++;
              S+=orderParameter(geo,p,lipids);
              A+=areaPerLipid(geo,p,lipids);
            }
          }

        /**
         * @brief Sample lipid order parameter
         *
         * @f[
         * S = \frac{1}{2} \left ( 3 (\mathbf{an})^2 -1 \right )
         * @\]
         *
         * where `a` is the unit vector between the tail and the head group,
         * `n` is the normal to the bilayer plane.
         */
        template<class Tcuboid, class Tpvec, class Tgroup>
          static double
          orderParameter(Tcuboid &geo, Tpvec &p, Tgroup &lipids, Point n=Point(0,0,1)) {
            Average<double> S;
            for (int i=0; i<lipids.numMolecules(); i++) {
              Group g;
              lipids.getMolecule(i,g); // i'th lipid
              Point a = geo.vdist( p[g.front()], p[g.back()]).normalized();
              S += 0.5 * ( 3 * pow(a.dot(n),2) - 1 );
            }
            return S.avg();
          }

        /**
         * @brief Sample area per lipid (normalized by sigma)
         */
        template<class Tcuboid, class Tpvec>
          static double
          areaPerLipid(Tcuboid &geo, Tpvec &p, Group &lipids) {
            return geo.len.x() * geo.len.y() / lipids.numMolecules() * 2
              / pow(2*p[ lipids.front() ].radius,2);
          }
    };

    /**
     * @brief Returns the dielectric constant outside the cutoff limit.
     *
     *        Only hold when using PBC and \f$\epsilon_{sur} = \epsilon\f$,
     *        [Neumann, M. (1983) Mol. Phys., 50, 841-858].
     *
     * @param pot The potential including geometry
     * @param spc The space including the particles
     * @param cutoff The cutoff of the reaction field
     */
    class getDielConst {
      private:
        Average<double> M;
        double volume;
        double convert;
        double cutoff;
      public:
        inline getDielConst(double cutoff_in) {
          cutoff = cutoff_in;
          convert = (3.33564*3.33564*(1e-30)/(0.20819434*0.20819434)); // Constant to convert to SI-units, including the cancelation of volume 10^-30
          volume = 4*pc::pi*pow(cutoff,3)/3;
          convert = convert*pc::pi/volume;
        }

        template<class Tpvec, class Tgeometry>
          void sample(const Tpvec &p, Tgeometry &geo) {
            Point origin(0,0,0), mu(0,0,0);
            for (auto &i : p)
              if (geo.sqdist(i,origin)<cutoff*cutoff)
                mu += i.mu*i.muscalar;
            M += mu.squaredNorm();
          }

        inline string info() {
          std::ostringstream o;
          if (M.cnt>0) {
            double Q = 0.25 + M.avg()*convert/pc::kT();
            o << "Eps: " << Q + std::sqrt(Q*Q+0.5) << "\n";
            //o << "<M>: " << M.avg() << ", convert/kT " << convert/pc::kT() << "\n";
          }
          return o.str();
        }
    };
    /*
    class checkWhiteNoise {
      private:
        std::vector<double> noise;
        int le;
        double significance;
        double lag;
        double mu;
        double sigma2;
      public:
        inline checkWhiteNoise(std::vector<double> noise_in,double significance_in, double lag_in) {
          noise = noise_in;
          significance = significance_in;
          lag = lag_in;
          le = noise.size();
          getMu();
          getVariance();
        }
        
        bool check(int h) {
          double Q = 0.0;
          for(int k = 0;k < h; k++) {
            Q += noise[k]*noise[k+lag]/(le-k);
          }
          Q = le*(le+2)*Q;
          double chi2 = getChi2(1-significance,h);
          if(Q > chi2) {
            return false;
          } else {
              return true;
            }
        }
        
        void getMu() {
          mu = 0;
          for(int k = 0;k < le; k++) {
             mu += noise[k];
          }
          mu /= le;
        }
        
        void getVariance() {
          sigma2 = 0;
          for(int k = 0;k < le; k++) {
             sigma2 += noise[k]*noise[k];
          }
          sigma2 -= mu*mu;
        }

        double getChi2(double x, int k) {
          if(x < 0) return 0.0;
          return (1-(incgamma(x,k)/(Gamma(x))));
        }
        
        double incgamma (double x, double a){
          double sum = 0;
          double term = 1.0/a;
          int n = 1;
          while (term != 0){
            sum = sum + term;
            term = term*(x/(a+n));
            n++;
          }
          return pow(x,a)*exp(-1*x)*sum;
        }

        double Gamma(double x) {
          if(std::abs(x-(int)x) < 1e-6 && (int)x > 1) {
              return (int)x*Gamma(((int)x)-1);
          } else if(std::abs(x-(int)x) < 1e-6 && (int)x == 1) {
            return 1;
          }
          
          if(x <= 0) { 
            return 0.0; 
          } else if(x <= 0.001) {
            double constant = 0.577215664901532860606512090; // Euler's gamma constant
            return 1.0/(x*(1.0 + constant*x));
          } else if(x <= 12) {
            double y = x;
            int n = 0;
            bool arg_was_less_than_one = (y < 1.0);
            if (arg_was_less_than_one) {
              y += 1.0;
            } else {
              n = static_cast<int> (floor(y)) - 1;  // will use n later
              y -= n;
            }
            static const double p[] =
            {
              -1.71618513886549492533811E+0,
              2.47656508055759199108314E+1,
              -3.79804256470945635097577E+2,
              6.29331155312818442661052E+2,
              8.66966202790413211295064E+2,
              -3.14512729688483675254357E+4,
              -3.61444134186911729807069E+4,
              6.64561438202405440627855E+4};
            static const double q[] =
            {
              -3.08402300119738975254353E+1,
              3.15350626979604161529144E+2,
              -1.01515636749021914166146E+3,
              -3.10777167157231109440444E+3,
              2.25381184209801510330112E+4,
              4.75584627752788110767815E+3,
              -1.34659959864969306392456E+5,
              -1.15132259675553483497211E+5};
            double num = 0.0;
            double den = 1.0;
            int i;
            double z = y - 1;
            for (i = 0; i < 8; i++) {
              num = (num + p[i])*z;
              den = den*z + q[i];
            }
            double result = num/den + 1.0;
            if (arg_was_less_than_one) {
              result /= (y-1.0);
            } else {
              for (i = 0; i < n; i++)
                result *= y++;
            }
            return result;
          } else if (x > 171.624)
            return pc::infty;
        }
    };*/
  }//namespace
}//namespace
#endif
