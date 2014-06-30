#ifndef FAUNUS_ANALYSIS_H
#define FAUNUS_ANALYSIS_H

#include <faunus/common.h>
#include <faunus/average.h>
#include <faunus/physconst.h>
#include <faunus/group.h>
#include <faunus/space.h>
#include <faunus/point.h>
#include <faunus/textio.h>
#include <faunus/energy.h>
#include <Eigen/Core>

#include <chrono>
#include <thread>

namespace Faunus {

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
     * @brief Pressure analysis using the virial theorem
     *
     * This calculates the excess pressure tensor defined as
     * @f[
     * \mathcal{P} = \frac{1}{3V}\left <
     * \sum_{i}^{N-1} \sum_{j=i+1}^N \mathbf{r}_{ij} \otimes \mathbf{f}_{ij}
     * \right >_{NVT}
     * @f]
     * where @f$r@f$ and @f$f@f$ are the distance and force, @f$V@f$ the system volume,
     * and the excess pressure scalar is the trace of @f$\mathcal{P}@f$.
     * The trivial kinetic contribution is currently not included.
     *
     * References:
     *
     * - <http://dx.doi.org/10/ffwrhd>
     * - <http://dx.doi.org/10/fspzcx>
     *
     * @todo At the moment this analysis is limited to "soft" systems only,
     * i.e. for non-rigid systems with continuous potentials.
     */
    class VirialPressure : public AnalysisBase {
      private:

        typedef Eigen::Matrix3d Ttensor;
        Ttensor T;           // excess pressure tensor
        Average<double> Pid; // ideal pressure

        inline string _info() {
          using namespace Faunus::textio;
          std::ostringstream o;

          if (cnt>0) {
            vector<double> P(3);
#ifdef __INTEL_COMPILER
            vector<string> id(3);
            id[0]="Ideal";
            id[1]="Excess";
            id[2]="Total";
#else
            vector<string> id = {"Ideal", "Excess", "Total"};
#endif
            P[0] = Pid.avg();       // ideal
            P[1] = (T/cnt).trace(); // excess
            P[2] = P[0] + P[1];     // total

            char l=15;
            double kT=pc::kB*pc::T();
            o << "\n  " << std::right
              << setw(l+l) << "kT/"+angstrom+cubed
              << setw(l) << "mM" << setw(l) << "Pa" << setw(l) << "atm" << "\n";
            for (int i=0; i<3; i++) {
              o << std::left << setw(l) << "  "+id[i] << std::right
                << setw(l) << P[i]
                << setw(l) << P[i] * 1e30/pc::Nav
                << setw(l) << P[i] * kT*1e30
                << setw(l) << P[i] * kT*1e30/0.980665e5
                << "\n";
            }
            o << "\n  Osmotic coefficient = " << 1+P[1]/P[0] << "\n";
            o << "  Excess pressure tensor (mM):\n\n"
              << T/cnt*1e30/pc::Nav << endl;
          }
          return o.str();
        }

        void _test(UnitTest &test) {
          test("virial_pressure_mM", (T/cnt).trace()*1e30/pc::Nav );
        }

        template<class Tpvec, class Tgeo, class Tpot>
          Ttensor g_internal(const Tpvec &p, Tgeo &geo, Tpot &pot, Group &g) { 
            Ttensor t;
            t.setZero();
            for (auto i=g.front(); i<g.back(); i++)
              for (auto j=i+1; j<=g.back(); j++) {
                auto rij = geo.vdist(p[i],p[j]);
                auto fij = pot.f_p2p(p[i],p[j]);
                t += rij * fij.transpose();
              }
            return t;
          }

        template<class Tpvec, class Tgeo, class Tpot>
          Ttensor g2g(const Tpvec &p, Tgeo &geo, Tpot &pot, Group &g1, Group &g2) { 
            assert(&g1!=&g2);
            Ttensor t;
            t.setZero();
            for (auto i : g1)
              for (auto j : g2) {
                auto rij = geo.vdist(p[i],p[j]);
                auto fij = pot.f_p2p(p[i],p[j]);
                t += rij * fij.transpose();
              }
            return t;
          }

      public:

        VirialPressure() {
          name="Virial Pressure";
          noMolecularPressure=false;
          T.setZero();
        }

        /*! @brief Ignore internal pressure in molecular groups (default: false) */
        bool noMolecularPressure;

        template<class Tspace, class Tpotential>
          void sample(Tspace &spc, Tpotential &pot, int d=3, double area=0) {
            cnt++;
            Ttensor t;
            t.setZero();

            int N=spc.p.size();
            double V=spc.geo.getVolume();
            if (d==2) {
              assert(area>0 && "Area must be specified for 2D sampling");
              V=area;
            }

            // loop over groups internally
            for (auto g : spc.groupList()) {
              if (noMolecularPressure)
                if (g->isMolecular()) {
                  N=N-g->size()+1;
                  continue;
                }
              t+=g_internal(spc.p, spc.geo, pot, *g);
            }

            // loop group-to-group
            auto beg=spc.groupList().begin();
            auto end=spc.groupList().end();
            for (auto gi=beg; gi!=end; ++gi)
              for (auto gj=gi; ++gj!=end;)
                t+=g2g(spc.p, spc.geo, pot, *(*gi), *(*gj));

            // add to grand avarage
            T += t/(d*V);
            Pid += N/V;
          }
    };

    /**
     * @brief General class for handling 2D tables - xy data, for example.
     * @date Lund 2011
     * @note `Tx` is used as the `std::map` key and which may be
     * problematic due to direct floating point comparison (== operator).
     * We have not experienced any issues with this, though. This uses
     * `std::map` and table lookup is of complexity logarithmic with N.
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
          template<class T=double>
            void save(string filename, T scale=1) {
              if (tabletype==HISTOGRAM) {
                if (!map.empty()) map.begin()->second*=2;   // compensate for half bin width
                if (map.size()>1) (--map.end())->second*=2; // -//-
              }

              if (!map.empty()) {
                std::ofstream f(filename.c_str());
                f.precision(10);
                if (f) {
                  for (auto m : map)
                    f << m.first << " " << get( m.first ) * scale << "\n";
                }
              }

              if (tabletype==HISTOGRAM) {
                if (!map.empty()) map.begin()->second/=2;   // restore half bin width
                if (map.size()>1) (--map.end())->second/=2; // -//-
              }
            }

          /** @brief Sums up all previous elements and saves table to disk */
          template<class T=double>
            void sumSave(string filename, T scale=1) {
              if (tabletype==HISTOGRAM) {
                if (!map.empty()) map.begin()->second*=2;   // compensate for half bin width
                if (map.size()>1) (--map.end())->second*=2; // -//-
              }

              if (!map.empty()) {
                std::ofstream f(filename.c_str());
                f.precision(10);
                if (f) {
                  Tx sum_t = 0.0;
                  for (auto m : map) {
                    sum_t += get( m.first );
                    f << m.first << " " << sum_t * scale << "\n";
                  }
                }
              }

              if (tabletype==HISTOGRAM) {
                if (!map.empty()) map.begin()->second/=2;   // restore half bin width
                if (map.size()>1) (--map.end())->second/=2; // -//-
              }
            }

          Tmap getMap() {
            return map;
          }

          Tx getResolution() {
            return dx;
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

          /**
           * @brief Convert table to matrix
           */
          Eigen::MatrixXd tableToMatrix() {
            assert(!this->map.empty() && "Map is empty!");
            Eigen::MatrixXd table(2,map.size());
            table.setZero();
            int I = 0;
            for (auto &m : this->map) {
              table(0,I) = m.first;
              table(1,I) = m.second;
              I++;
            }
            return table;
          }
      };

    /**
     * @brief Subtract two tables
     */
    template<class Tx, class Ty, class Tmap>
      Table2D<Tx,Ty> operator-(Table2D<Tx,Ty> &a, Table2D<Tx,Ty> &b) {
        assert(a.tabletype=b.tabletype && "Table a and b needs to be of same type");
        Table2D<Tx,Ty> c(std::min(a.getResolution(),b.getResolution()),a.tabletype);
        Tmap a_map = a.getMap();
        Tmap b_map = b.getMap();

        if (a.tabletype=="HISTOGRAM") {
          if (!a_map.empty()) a_map.begin()->second*=2;   // compensate for half bin width
          if (a_map.size()>1) (--a_map.end())->second*=2; // -//-
          if (!b_map.empty()) b_map.begin()->second*=2;   // compensate for half bin width
          if (b_map.size()>1) (--b_map.end())->second*=2; // -//-
        }

        for (auto &m1 : a_map) {
          for (auto &m2 : b_map) {
            c(m1.first) = m1.second-m2.second;
            break;
          }
        }

        if (a.tabletype=="HISTOGRAM") {
          if (!a_map.empty()) a_map.begin()->second/=2;   // compensate for half bin width
          if (a_map.size()>1) (--a_map.end())->second/=2; // -//-
          if (!b_map.empty()) b_map.begin()->second/=2;   // compensate for half bin width
          if (b_map.size()>1) (--b_map.end())->second/=2; // -//-
        }
        return c;
      }

    /**
     * @brief Addition two tables
     */
    template<class Tx, class Ty, class Tmap>
      Table2D<Tx,Ty> operator+(Table2D<Tx,Ty> &a, Table2D<Tx,Ty> &b) {
        assert(a.tabletype=b.tabletype && "Table a and b needs to be of same type");
        Table2D<Tx,Ty> c(std::min(a.getResolution(),b.getResolution()),a.tabletype);
        Tmap a_map = a.getMap();
        Tmap b_map = b.getMap();

        if (a.tabletype=="HISTOGRAM") {
          if (!a_map.empty()) a_map.begin()->second*=2;   // compensate for half bin width
          if (a_map.size()>1) (--a_map.end())->second*=2; // -//-
          if (!b_map.empty()) b_map.begin()->second*=2;   // compensate for half bin width
          if (b_map.size()>1) (--b_map.end())->second*=2; // -//-
        }

        for (auto &m : a_map) {
          c(m.first) += m.second;
        }
        for (auto &m : b_map) {
          c(m.first) += m.second;
        }

        if (a.tabletype=="HISTOGRAM") {
          if (!a_map.empty()) a_map.begin()->second/=2;   // compensate for half bin width
          if (a_map.size()>1) (--a_map.end())->second/=2; // -//-
          if (!b_map.empty()) b_map.begin()->second/=2;   // compensate for half bin width
          if (b_map.size()>1) (--b_map.end())->second/=2; // -//-
        }

        return c;
      }

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
            return 4./3.*pc::pi*( pow(x+0.5*this->dx,3)
                - pow(x-0.5*this->dx,3) );
          }

          double get(Tx x) {
            assert( volume(x)>0 );
            assert( this->count()>0 );

            if (bulkconc.cnt==0) bulkconc+=1;
            if (bulkconc.avg()<1e-6) bulkconc+=1;
            if (Npart.cnt==0) Npart+=1;

            return ((double)this->operator()(x)*Npart.avg()) / (volume(x) *(double)this->count() * bulkconc.avg())
              ;

          }
          Average<double> bulkconc; //!< Average bulk concentration
          Average<double> Npart;
        public:
          Tx maxdist; //!< Pairs with distances above this value will be skipped (default: infinity)

          /*!
           * \param res Resolution of X axis
           */
          RadialDistribution(Tx res=0.2) : Ttable(res,Ttable::HISTOGRAM) {
            this->name="Radial Distribution Function";

            maxdist=pc::infty;
            static_assert( std::is_integral<Ty>::value,
                "Histogram must be of integral type");
            static_assert( std::is_unsigned<Ty>::value,
                "Histogram must be unsigned");
          }
          /*!
           * \brief Sample radial distibution of two atom types
           * \param spc Simulation space
           * \param g Group to search
           * \param ida Atom id of first particle
           * \param idb Atom id of second particle
           */
          template<class Tspace>
            void sample(Tspace &spc, Group &g, short ida, short idb) {
              for (auto i=g.begin(); i!=g.end()-1; i++)
                for (auto j=i+1; j!=g.end(); j++)
                  if ( (spc.p[*i].id==ida && spc.p[*j].id==idb) || (spc.p[*i].id==idb && spc.p[*j].id==ida) ) {
                    Tx r=spc.geo.dist(spc.p[*i], spc.p[*j]);
                    if (r<=maxdist)
                      this->operator() (r)++;
                  }
              int bulk=0;
              for (auto i : g){
                if (spc.p[i].id==ida || spc.p[i].id==idb){
                  bulk++;
                }
              }
              Npart+=bulk;
              bulkconc += bulk / spc.geo.getVolume();
            }


          template<class Tspace>
            void sample(Tspace &spc, short ida, short idb) {
              Group all(0, spc.p.size()-1);
              assert(all.size()==(int)spc.p.size());
              return sample(spc,all,ida,idb);
            }

          template<class Tspace>
            void sampleMolecule(Tspace &spc, Group &sol) {
              for (int i=0; i<sol.numMolecules()-1; i++) {
                for (int j=i+1; j<sol.numMolecules(); j++) {
                  Group ig, jg;
                  sol.getMolecule(i,ig);
                  sol.getMolecule(j,jg);
                  Point icm = ig.massCenter(spc);
                  Point jcm = jg.massCenter(spc);
                  this->operator() (spc.geo.dist(icm,jcm))++;
                }
              }
            }

          // Same as sampeMolecule but different inputs
          template<class Tspace>
            void sampleMoleculeGroup(Tspace &spc, vector<Group> &g, string name) {
              int bulk = 0;
              for(size_t i = 0; i < g.size()-1; i++) {
                Group ig = g[i];
                if(ig.name == name) {
                  bulk++;
                  for(size_t j = i+1; j < g.size(); j++) {
                    Group jg = g[j];
                    if(jg.name == name) {
                      Point icm = ig.massCenter(spc);
                      Point jcm = jg.massCenter(spc);
                      this->operator() (spc.geo.dist(icm,jcm))++;
                    }
                  }
                }
              }
              Npart+=bulk;
              bulkconc += bulk / spc.geo.getVolume();
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

    /**
     * @brief Base class for force calculations
     *
     * Includes some neccessary functionality for deriving the force.
     *
     * @author Axel Thuresson
     * @date Lund, 2013
     */
    class TwobodyForce : public AnalysisBase {
      protected:
        string _info();         //!< Print results of analysis
        Group* igroup1;
        Group* igroup2;
        Group* ions;
      public:
        virtual Point meanforce();
        TwobodyForce(InputMap&, Group&, Group&, Group&);//!< Constructor
        void save(string);
        void setTwobodies(Group&, Group&, Group&);
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
        TwobodyForceDirect(InputMap&, Group&, Group&, Group&);//!< Constructor
        Point meanforce();

        /** @brief Calculate the direct force between the two bodies */
        template<class Tpvec, class Tenergy>
          void calc(Tpvec &p, Tenergy &pot) {
            if (run()) {
              // Force between the two bodies
              for (auto i : *igroup1) {
                for (auto j : *igroup2) {
                  Point f = pot.f_p2p(p[i], p[j]);
                  f_pp += f;
                }
              }
              //f_pp += 1.0*_f_pp;
              //f_mean1 += 1.0*_f_pp;
              //f_mean2 += -1.0*_f_pp;
              for (auto i : *igroup1) {
                for (auto j : *ions) {
                  Point f = pot.f_p2p(p[i],p[j]);
                  f_pi += f;
                }
              }
              for (auto i : *igroup2) {
                for (auto j : *ions) {
                  Point f = pot.f_p2p(p[i], p[j]);
                  f_ip += f;
                }
              }
            }
          }
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
        TwobodyForceMidp(InputMap&, Group&, Group&, Group&, Analysis::LineDistributionNorm<float,unsigned long int>*);//!< Constructor
        Point meanforce();

        /** @brief Calculate the direct force between the two bodies */
        template<class Tpvec, class Tenergy>
          void calc(Tpvec &p, Tenergy &pot) {
            if (run()) {
              // Force between the two bodies
              for (auto i : *igroup1) {
                for (auto j : *igroup2) {
                  Point f = pot.f_p2p(p[i], p[j]);
                  f_pp += f;
                }
              }

              for (auto i : *igroup1) {
                for (auto j : *ions) {
                  if (p[j].z() < 0.0) {
                    Point f = pot.f_p2p(p[i], p[j]);
                    f_pi += f;
                  }
                }
              }

              for (auto i : *igroup2) {
                for (auto j : *ions) {
                  if (p[j].z() >= 0.0) {
                    Point f = pot.f_p2p(p[i], p[j]);
                    f_ip += f;
                  }
                }
              }

              for (auto i : *ions) {
                if (p[i].z() >= 0.0) {
                  for (auto j : *ions) {
                    if (p[j].z() < 0.0) {
                      Point f = pot.f_p2p(p[i],p[j]);
                      f_ii += f;
                    }
                  }
                }
              }
            }
          }
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
            assert( spc.geo.dist(pol.cm, pol.massCenter(spc))<1e-9
                && "Mass center must be in sync.");
            Point rg2=vectorgyrationRadiusSquared(pol,spc);
            return rg2.x()+rg2.y()+rg2.z();
          }

        template<class Tgroup, class Tspace>
          Point vectorEnd2end(const Tgroup &pol, const Tspace &spc) {
            return spc.geo.vdist( spc.p[pol.front()], spc.p[pol.back()] );
          }

      public:
        PolymerShape();

        /** 
         * @note This functions is now public and const. I don't see the point of making it static, yet. - Joao Henriques.
         */
        template<class Tgroup, class Tspace>
          Point vectorgyrationRadiusSquared(const Tgroup &pol, const Tspace &spc) const {
            assert( spc.geo.dist(pol.cm, pol.massCenter(spc))<1e-9
                && "Mass center must be in sync.");
            double sum=0;
            Point t, r2(0,0,0);
            for (auto i : pol) {
              t = spc.p[i]-pol.cm;                     // vector to center of mass
              spc.geo.boundary(t);                     // periodic boundary (if any)
              r2.x() += spc.p[i].mw * t.x() * t.x();
              r2.y() += spc.p[i].mw * t.y() * t.y();
              r2.z() += spc.p[i].mw * t.z() * t.z();
              sum += spc.p[i].mw;                      // total mass
            }
            assert(sum>0 && "Zero molecular weight not allowed.");
            return r2*(1./sum);
          }

        /** @brief Sample properties of Group (identified by group name) */
        template<class Tgroup, class Tspace>
          void sample(const Tgroup &pol, const Tspace &spc) {
            if (!run() || pol.front()==pol.back())
              return;
            Point r2 = vectorgyrationRadiusSquared(pol,spc);
            double rg2 = r2.x()+r2.y()+r2.z(); 
            double re2 = spc.geo.sqdist( spc.p[pol.front()], spc.p[pol.back()] );
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
            assert( spc.geo.dist(g.cm, g.massCenter(spc))<1e-9
                && "Mass center must be in sync.");
            Point t, mu(0,0,0);
            for (auto i : g)
              if (!excluded(spc.p[i])) {
                t = spc.p[i]-g.cm;                // vector to center of mass
                spc.geo.boundary(t);               // periodic boundary (if any)
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
     * @brief Multipolar decomposition between groups as a function of separation
     *
     * This will analyse the electrostatic energy between two groups as
     * a function of their mass center separation. Sampling consists of
     * the following:
     *
     * 1. The exact electrostatic energy is calculated by explicitly summing
     *    Coulomb interactions between charged particles
     * 2. Each group -- assumed to be a molecule -- is translated into a
     *    multipole (monopole, dipole, quadrupole)
     * 3. Multipolar interaction energies are calculated, summed, and tabulated
     *    together with the exact electrostatic interaction energy. Ideally
     *    (infinite number of terms) the multipoles should capture full
     *    electrostatics
     *
     * The points 1-3 above will be done as a function of group-to-group
     * mass center separation.
     *
     * Note also that the moments are defined with
     * respect to the *mass* center, not *charge* center. While for most
     * macromolecules there is only a minor difference between the two,
     * the latter is more appropriate and is planned for a future update.
     * A simply way to mimic this is to assign zero mass to all neutral
     * atoms in the molecule.
     *
     * The constructor takes the following `InputMap` keywords:
     *
     * Keyword             | Note
     * :------------------ | :--------------------------------------------
     * `multipoledist_dr`  | Distance resolution - default is 0.2 angstrom
     *
     * @date Malmo 2014
     * @note Needs testing!
     * @todo Add option to use charge center instead of mass center
     */
    template<class Tspace, class Tcoulomb=Potential::Coulomb>
      class MultipoleDistribution : public AnalysisBase {
        private:
          Energy::Nonbonded<Tspace,Tcoulomb> pot; // Coulombic hamiltonian

          string _info() FOVERRIDE {
            return textio::header("Multipole Analysis");
          }

          struct data {
            double tot, ii, id, iq, dd;
            unsigned long int cnt;
            data() : cnt(0) {}
          };

          std::map<int,data> m; // slow, but OK for infrequent sampling
          double dr;            // distance resolution

          template<class Tgroup>
            Tensor<double> quadrupoleMoment(const Tspace &s, const Tgroup &g) const {
              Tensor<double> theta;
              theta.clear();
              assert(g.size()<=(int)s.p.size());
              for (auto i : g) {
                Point t=s.p[i] - g.cm;
                s.geo.boundary(t);
                theta = theta + t*t.transpose()*s.p[i].charge;
              }
              return 0.5*theta;
            }

          // convert group to multipolar particle
          template<class Tgroup>
            DipoleParticle toMultipole(const Tspace &spc, const Tgroup &g) const {
              DipoleParticle m;
              m = g.cm;
              m.charge = netCharge(spc.p, g);            // monopole
              m.mu = Geometry::dipoleMoment(spc, g); // dipole
              m.muscalar = m.mu.norm();
              if (m.muscalar>1e-8)
                m.mu=m.mu/m.muscalar; 
              m.theta = quadrupoleMoment(spc, g);    // quadrupole
              return m;
            }

          template<class Tgroup>
            double g2g(const Tspace &spc, const Tgroup &g1, const Tgroup &g2) { 
              double u=0;
              for (auto i : g1)
                for (auto j : g2)
                  u += spc.p[i].charge * spc.p[j].charge /
                    spc.geo.dist( spc.p[i], spc.p[j] );
              return u * pot.pairpot.bjerrumLength();
            }

        public:

          MultipoleDistribution(InputMap &in) : pot(in) {
            dr = in.get("multipoledist_dr", 0.1,
                "Distance resolution of multipole analysis");
          }

          /**
           * @brief Sample multipole energy
           * @param spc Simulation space
           * @param g1  Group with molecule 1
           * @param g2  Group with molecule 2
           * @note Group mass-centers (`Group::cm`) must be up-to-date before
           *       calling this function
           */
          template<class Tmultipole=DipoleParticle>
            void sample(Tspace &spc, Group &g1, Group &g2) {
              if (run()) {
                // multipoles and cm-cm distance
                auto a = toMultipole(spc,g1);
                auto b = toMultipole(spc,g2);
                auto r = spc.geo.vdist(g1.cm, g2.cm);
                double r2inv = 1/r.squaredNorm();
                double rinv = sqrt(r2inv);
                double r3inv = rinv * r2inv;

                // multipolar energy
                pot.setSpace(spc);
                data d;
                d.cnt++;
                d.tot = pot.g2g(spc.p, g1, g2); // exact el. energy
                d.ii = a.charge * b.charge * rinv; // ion-ion, etc.
                d.id = ( a.charge*b.mu.dot(r) - b.charge*a.mu.dot(r) ) * r3inv;
                d.dd = mu2mu(a.mu, b.mu, a.muscalar*b.muscalar, r);
                d.iq = q2quad(a.charge, b.theta, b.charge, a.theta, r);

                // add to grand average
                int key = to_bin(1/rinv,dr);
                auto it = m.find(key);
                if (it==m.end())
                  m[key]=d;
                else {
                  it->second.cnt++;
                  it->second.ii  += d.ii;
                  it->second.id  += d.id;
                  it->second.iq  += d.iq;
                  it->second.dd  += d.dd;
                  it->second.tot += d.tot;
                }
              }
            }

          /** @brief Save multipole distribution to disk */
          void save(const string &filename) const {
            std::ofstream f(filename.c_str());
            if (f) {
              char w=12;
              auto lB=pot.pairpot.bjerrumLength();
              f.precision(4);
              f << "# Multipolar energy analysis (kT)\n"
                << std::left << setw(w) << "# r/AA" << std::right << setw(w) << "exact"
                << setw(w) << "total" << setw(w) << "ionion" << setw(w) << "iondip"
                << setw(w) << "dipdip" << setw(w) << "ionquad\n";
              for (auto &i : m)
                f << std::left
                  << setw(w) << i.first*dr                  // r
                  << std::right
                  << setw(w) << i.second.tot/i.second.cnt   // exact (already in kT)
                  << setw(w) << lB*(i.second.ii+i.second.id+i.second.dd+i.second.iq)/i.second.cnt // total
                  << setw(w) << lB*i.second.ii/i.second.cnt // individual poles...
                  << setw(w) << lB*i.second.id/i.second.cnt
                  << setw(w) << lB*i.second.dd/i.second.cnt
                  << setw(w) << lB*i.second.iq/i.second.cnt
                  << "\n";
            }
          }
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

          void add(Tparticle p) { g.push_back(p); }

          /* @brief Add particle to insert - sum of added particle charges should be zero.*/
          template<class Tpvec>
            void add(Tpvec &p) {
              std::map<short,bool> map; // replace w. `std::set`
              for (auto i : p)
                map[ i.id ] = true;
              for (auto &m : map) {
                Tparticle a;
                a=atom[m.first];
                add(a);
              }
            }

          /** @brief Sampled mean activity coefficient */
          double gamma() { return exp(muex()); }

          /** @brief Sampled mean excess chemical potential */
          double muex() { return -log(expsum.avg())/g.size() ; }

          /** @brief Insert and analyse `n` times */
          template<class Tspace, class Tenergy>
            void sample(Tspace &spc, Tenergy &pot, int ghostin) {
              int n=g.size();
              if (n>0)
                if (run()) {
                  while (ghostin-->0) {
                    double du=0;
                    for (auto &i : g)
                      spc.geo.randompos(i);     // random ghost positions
                    for (auto &i : g)
                      du+=pot.all2p(spc.p, i);  // energy with all particles in space
                    for (int i=0; i<n-1; i++)
                      for (int j=i+1; j<n; j++)
                        du+=pot.p2p(g[i], g[j]);// energy between ghost particles
                    expsum += exp(-du);
                  }
                }
            }
      };

    /**
     * @brief Single particle hard sphere Widom insertion with charge scaling
     *
     * This will calculate excess chemical potentials for single particles
     * in the primitive model of electrolytes. Use the `add()` functions
     * to add test or *ghost* particles and call `sample()` to perform single
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
     * @warning Works only for the primitive model
     * @note This is a conversion of the Widom routine found in the `bulk.f`
     *       fortran program by Bolhuis/Jonsson/Akesson at Lund University.
     * @author Martin Trulsson and Mikael Lund
     * @date Lund / Prague 2007-2008.
     */
    template<class Tparticle>
      class WidomScaled : public AnalysisBase {

        private:

          typedef std::vector<double> Tvec;
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
          int ghostin;        //< ghost insertions
          double lB;          //!< Bjerrum length [a]

          void init() {
            int gspec=g.size();
            chel.resize(gspec);
            chhc.resize(gspec);
            chex.resize(gspec);
            chtot.resize(gspec);
            ewden.resize(gspec);
            ewnom.resize(gspec);
            chint.resize(gspec);
            expuw.resize(gspec);
            chexw.resize(gspec);
            ihc.resize(gspec);
            irej.resize(gspec);

            for (int i=0; i<gspec; i++){
              chel[i]=0;
              chhc[i]=0;
              chex[i]=0;
              chtot[i]=0;
              ihc[i]=0;
              ewden[i].resize(11);
              ewnom[i].resize(11);
              chint[i].resize(11);
              for(int j=0; j<11; j++)
                ewden[i][j] = ewnom[i][j] = chint[i][j] = 0;
            }
          }

          template<class Tgeo>
            bool overlap(const Tparticle &a, const Tparticle &b, const Tgeo &geo)
            {
              double s=a.radius+b.radius;
              return (geo.sqdist(a,b)<s*s) ? true : false;
            }

          string _info() {
            using namespace textio;
            std::ostringstream o;
            double aint4,aint2,aint1;
            for(size_t i=0; i<g.size(); i++) {
              for(int cint=0; cint<11; cint++) {
                if(ewden[i][cint]==0)
                  std::cerr << "# WARNING: Widom denominator equals zero" << endl;
                else
                  chint[i][cint]=ewnom[i][cint]/ewden[i][cint];
              }
              aint4=chint[i][1]+chint[i][3]+chint[i][5]+chint[i][7]+chint[i][9];
              aint2=chint[i][2]+chint[i][4]+chint[i][6]+chint[i][8];
              aint1=chint[i][0]+chint[i][10];  
              chel[i]=1./30.*(aint1+2*aint2+4*aint4);
            }

            int cnttot=cnt*ghostin;
            o << pad(SUB,w,"Number of insertions") << cnttot << endl
              << pad(SUB,w,"Excess chemical potentials (kT)") << endl
              << "             total    elec.   hs             z        r"<< endl;
            char w=10;
            for (size_t i=0; i < g.size(); i++) {
              chhc[i]=-log(double(cnttot-ihc[i])/cnttot);
              chexw[i]=-log(expuw[i]);
              chex[i]=chhc[i]+chel[i];
              o.unsetf( std::ios_base::floatfield );
              o << "    [" << i << "] "
                << std::setprecision(4)
                << std::setw(w) << chex[i]
                << std::setw(w) << chel[i]
                << std::setw(w) << chhc[i]
                << std::setprecision(2) << std::fixed
                << std::setw(w) << g[i].charge
                << std::setw(w) << g[i].radius << endl;
            }
            return o.str();
          }

        public:

          /**
           * @param bjerrumLength Bjerrum length [angstrom]
           * @param insertions Number of insertions per call to `insert()`
           */
          WidomScaled(double bjerrumLength, int insertions) {
            assert(insertions>=0);
            assert(bjerrumLength>0);
            name="Single particle Widom insertion w. charge scaling"; 
            cite="doi:10/ft9bv9 + doi:10/dkv4s6"; 
            lB=bjerrumLength;
            ghostin=insertions;
          }

          /**
           * @brief Add ghost particle
           *
           * This will add particle `p` to the list of ghost particles
           * to insert.
           */
          void add(const Tparticle &p) {
            g.push_back(p);
            init();
          }

          /**
           * @brief Add ghost particles
           *
           * This will scan the particle vector for particles and each unique type
           * will be added to the list a ghost particles to insert.
           */
          template<class Tpvec>
            void add(const Tpvec &p) {
              std::set<particle::Tid> ids;
              for (auto &i : p)
                ids.insert(i.id);
              for (auto i : ids) {
                particle a;
                a=atom[i];
                add(a);
              }
            }

          /**
           * @brief Do test insertions and sample excess chemical potential 
           *
           * @param p List of particles to insert into. This will typically be the main
           *          particle vector, i.e. `Space::p`.
           * @param geo Geometry to use for distance calculations and random position generation
           */
          template<class Tpvec, class Tgeo>
            void sample(const Tpvec &p, Tgeo &geo) {
              assert(lB>0);
              assert(&geo!=nullptr);
              if (!g.empty())
                if (!p.empty())
                  if (run()) {
                    Tparticle ghost;
                    double u,cu;
                    for (int i=0; i<ghostin; i++) {
                      geo.randompos(ghost);
                      int goverlap=0;
                      for (size_t k=0; k<g.size(); k++) {
                        ghost.radius = g[k].radius;
                        irej[k]=0;
                        int j=0;
                        while (!overlap(ghost,p[j],geo) && j<(int)p.size())
                          j++;
                        if (j!=(int)p.size()) {
                          ihc[k]++;
                          irej[k]=1;
                          goverlap++;
                        }
                      }

                      if ( goverlap != (int)g.size() ) {
                        cu=0;
                        u=0;  //elelectric potential (Coulomb only!)
                        for (auto &i : p) {
                          double invdi=1/geo.dist(ghost,i);
                          cu+=invdi;
                          u+=invdi*i.charge;
                        } 
                        cu=cu*lB;
                        u=u*lB;
                        double ew,ewla,ewd;
                        for (size_t k=0; k < g.size(); k++) {
                          if (irej[k]==0) {
                            expuw[k]+=exp(-u*g[k].charge);
                            for (int cint=0; cint<11; cint++) {
                              ew=g[k].charge*(u-double(cint)*0.1*g[k].charge*cu/double(p.size()));
                              ewla = ew*double(cint)*0.1;
                              ewd=exp(-ewla);
                              ewden[k][cint]+=ewd;
                              ewnom[k][cint]+=ew*ewd;
                            }
                          }
                        }
                      }
                    }
                  } 
            }

      }; // end of WidomScaled

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
     * @brief Analyse dielectric constant outside the cutoff limit.
     *
     * @note [Neumann, M. (1983) Mol. Phys., 50, 841-858].
     *
     * @param spc The space
     * @param filename Extention of filename from previous saved run (optional)
     */
    class DipoleAnalysis {
      private:
        Analysis::RadialDistribution<> rdf;
        Analysis::Table2D<double,double> kw, mucorr_angle;
        Analysis::Table2D<double,Average<double> > mucorr, mucorr_dist; 
        Analysis::Histogram<double,unsigned int> HM_x,HM_y,HM_z,HM_x_box,HM_y_box,HM_z_box,HM2,HM2_box;
        Average<double> M_x,M_y,M_z,M_x_box,M_y_box,M_z_box,M2,M2_box,diel_std, V_t;
        Average<double> mu_abs;
        
        int sampleKW;
        double cutoff2, N;
        double const_Diel, const_DielTinfoil, const_DielCM, const_DielRF, epsRF, constEpsRF;

      public:
        template<class Tspace, class Tinputmap>
          DipoleAnalysis(const Tspace &spc, Tinputmap &in) : rdf(0.1),kw(0.1),mucorr_angle(0.1),mucorr(0.1),mucorr_dist(0.1),HM_x(0.1),HM_y(0.1),HM_z(0.1),HM_x_box(0.1),HM_y_box(0.1),HM_z_box(0.1),HM2(0.1),HM2_box(0.1) {
            sampleKW = 0;
            setCutoff(spc.geo.len_half.x());
            N = spc.p.size();
            const_Diel = pc::e*pc::e*1e10/(3*pc::kT()*pc::e0);
            updateConstants(spc.geo.getVolume());
            updateDielectricConstantRF(in.get("epsilon_rf",80.));
            load(in.get("dipole_data_ext", string("")));
          }

        void setCutoff(double cutoff) {
          cutoff2 = cutoff*cutoff;
        }
        
        void updateDielectricConstantRF(double er) {
          epsRF = er;
          constEpsRF = 2*(epsRF - 1)/(2*epsRF + 1);
        }
        
        void updateConstants(double volume) {
          V_t += volume;
          const_DielTinfoil = const_Diel/V_t.avg();
          const_DielCM = const_DielTinfoil/3;
          const_DielRF = const_DielTinfoil/3;
        }

        /**
         * @brief Samples dipole-moment from dipole particles.
         * @param spc The space.
         */
        template<class Tspace>
          void sampleDP(const Tspace &spc) {
            Point origin(0,0,0);
            Point mu(0,0,0);        // In e\AA
            Point mu_box(0,0,0);    // In e\AA

            for (auto &i : spc.p) {
              if (spc.geo.sqdist(i,origin)<cutoff2) {
                mu += i.mu*i.muscalar;
              } else {
                mu_box += i.mu*i.muscalar;
              }
              mu_abs += i.muscalar;
            }
            mu_box += mu;
            samplePP(spc,mu,mu_box);
          }

        /**
         * @brief Samples dipole-moment from point particles.
         * 
         * @param spc The space
         * @param mu Dipoles to add to from within cutoff (optional)
         * @param mu_box Dipoles to add to from entire box (optional)
         */
        template<class Tspace>
          void samplePP(Tspace &spc, Point mu=Point(0,0,0), Point mu_box=Point(0,0,0)) {
            updateConstants(spc.geo.getVolume());
            Group all(0,spc.p.size()-1);
            all.setMassCenter(spc);
            mu += Geometry::dipoleMoment(spc,all,sqrt(cutoff2));
            mu_box += Geometry::dipoleMoment(spc,all);

            HM_x(mu.x())++;
            HM_y(mu.y())++;
            HM_z(mu.z())++;
            HM_x_box(mu_box.x())++;
            HM_y_box(mu_box.y())++;
            HM_z_box(mu_box.z())++;
            M_x += mu.x();
            M_y += mu.y();
            M_z += mu.z();
            M_x_box += mu_box.x();
            M_y_box += mu_box.y();
            M_z_box += mu_box.z();
            double sca = mu.dot(mu);
            HM2(sca)++;
            M2 += sca;
            sca = mu_box.dot(mu_box);
            HM2_box(sca)++;
            M2_box += sca;
            diel_std.add(getDielTinfoil());
          }

        /**
         * @brief Samples g(r), \f$ <\mu(0) \cdot \mu(r)> \f$, \f$ <\frac{1}{2} ( 3 \mu(0) \cdot \mu(r) - 1 )> \f$, Histogram(\f$ \mu(0) \cdot \mu(r) \f$) and distant-dependent Kirkwood-factor.
         * 
         * @param spc The space
         *
         */
        template<class Tspace>
          void sampleMuCorrelationAndKirkwood(Tspace &spc) {
            double r, sca;
            int N = spc.p.size() - 1;
            for (int i = 0; i < N; i++) {
              kw(0) += spc.p[i].mu.dot(spc.p[i].mu)*spc.p[i].muscalar*spc.p[i].muscalar;
              for (int j = i+1.; j < N + 1; j++) {
                r = spc.geo.dist(spc.p[i],spc.p[j]); 
                rdf(r)++;
                sca = spc.p[i].mu.dot(spc.p[j].mu);
                mucorr_angle(sca) += 1.;
                mucorr(r) += sca;
                mucorr_dist(r) += 0.5*(3*sca*sca -1.);
                kw(r) += 2*sca*spc.p[i].muscalar*spc.p[j].muscalar;
              }
            }
            kw(0) += spc.p[N].mu.dot(spc.p[N].mu)*spc.p[N].muscalar*spc.p[N].muscalar;
            sampleKW++;
          }

        /**
         * @brief Returns dielectric constant (using Tinfoil conditions).
         * \f$ 1 + \frac{<M^2>}{3V\epsilon_0k_BT} \f$
         */
        double getDielTinfoil() {
          // 1 + ( ( ( 4 * pi * <M^2> ) / ( 3 * V * kT ) ) / ( 4 * pi * e0 ) )
          return ( 1 + M2_box.avg()*const_DielTinfoil); 
        }  

        /**
         * @brief Returns dielectric constant ( Clausius-Mossotti ).
         * \f$ \frac{1 + \frac{2<M^2>}{9V\epsilon_0k_BT}}{1 - \frac{<M^2>}{9V\epsilon_0k_BT}} \f$
         * 
         * @warning Needs to be checked!
         */
        double getDielCM() {
          double temp = M2_box.avg()*const_DielCM;
          return ((1 + 2*temp)/(1 - temp));
        }  
        
        /**
         * @brief Returns dielectric constant ( Reaction Field ).
         * 
         * @warning Needs to be checked!
         */
        double getDielRF() {
          double avgRF = M2_box.avg()*const_DielRF;
          double sqrtRF = 0;
          if(-4*constEpsRF*avgRF + 1 > 0)
            sqrtRF = 3*sqrt(-4*constEpsRF*avgRF + 1);
          double one = 0.5*(2*constEpsRF - 4*avgRF + 1 + sqrtRF)/(constEpsRF + avgRF - 1);
          //double two = 0.5*(2*constEpsRF - 4*avgRF + 1 - sqrtRF)/(constEpsRF + avgRF - 1);
          return one;
        }  
        
        /**
         * @brief Saves data to files. 
         * @param nbr Extention of filename
         * 
         * @note \f$ g(r) \rightarrow \f$ gofr.dat+nbr
         *       \f$ \mu(0)\cdot\mu(r) \rightarrow \f$ mucorr.dat+nbr
         *       \f$ <\frac{1}{2} ( 3 \mu(0) \cdot \mu(r) - 1 )> \rightarrow \f$ mucorr_dist.dat+nbr
         * 
         */
        void save(string ext="") {
          rdf.save("gofr.dat"+ext);
          mucorr.save("mucorr.dat"+ext);
          mucorr_angle.save("mucorr_angle.dat"+ext);
          mucorr_dist.save("mucorr_dist.dat"+ext);
          kw.sumSave("kirkwood.dat"+ext,1.0/double(sampleKW));
          HM_x.save("hist_dip_x.dat"+ext);
          HM_y.save("hist_dip_y.dat"+ext);
          HM_z.save("hist_dip_z.dat"+ext);
          HM_x_box.save("hist_dip_x_box.dat"+ext);
          HM_y_box.save("hist_dip_y_box.dat"+ext);
          HM_z_box.save("hist_dip_z_box.dat"+ext);
          HM2.save("hist_dip2.dat"+ext);
          HM2_box.save("hist_dip2_box.dat"+ext);

          string filename = "dipoledata.dat"+ext;
          std::ofstream f(filename.c_str());
          f.precision(10);
          if (f) {
            if (M_x.cnt != 0)      f << "M_x " << M_x.cnt << " " << M_x.avg() << " " << M_x.sqsum << "\n";
            if (M_y.cnt != 0)      f << "M_y " << M_y.cnt << " " << M_y.avg() << " " << M_y.sqsum << "\n";
            if (M_z.cnt != 0)      f << "M_z " << M_z.cnt << " " << M_z.avg() << " " << M_z.sqsum << "\n";
            if (M_x_box.cnt != 0)  f << "M_x_box " << M_x_box.cnt << " " << M_x_box.avg() << " " << M_x_box.sqsum << "\n";
            if (M_y_box.cnt != 0)  f << "M_y_box " << M_y_box.cnt << " " << M_y_box.avg() << " " << M_y_box.sqsum << "\n";
            if (M_z_box.cnt != 0)  f << "M_z_box " << M_z_box.cnt << " " << M_z_box.avg() << " " << M_z_box.sqsum << "\n";
            if (M2.cnt != 0)       f << "M2 " << M2.cnt << " " << M2.avg() << " " << M2.sqsum << "\n";
            if (M2_box.cnt != 0)   f << "M2_box " << M2_box.cnt << " " << M2_box.avg() << " " << M2_box.sqsum << "\n";
            if (diel_std.cnt != 0) f << "diel_std " << diel_std.cnt << " " << diel_std.avg() << " " << diel_std.sqsum;
            if (mu_abs.cnt != 0)      f << "mu_abs " << mu_abs.cnt << " " << mu_abs.avg() << " " << mu_abs.sqsum << "\n";
          }
        }

        void load(string ext="") {
          if(ext == "none")
            ext = "";
          rdf.load("gofr.dat"+ext);
          mucorr.load("mucorr.dat"+ext);
          mucorr_angle.load("mucorr_angle.dat"+ext);
          mucorr_dist.load("mucorr_dist.dat"+ext);
          kw.load("kirkwood.dat"+ext);
          HM_x.load("hist_dip_x.dat"+ext);
          HM_y.load("hist_dip_y.dat"+ext);
          HM_z.load("hist_dip_z.dat"+ext);
          HM_x_box.load("hist_dip_x_box.dat"+ext);
          HM_y_box.load("hist_dip_y_box.dat"+ext);
          HM_z_box.load("hist_dip_z_box.dat"+ext);
          HM2.load("hist_dip2.dat"+ext);
          HM2_box.load("hist_dip2_box.dat"+ext);

          string filename = "dipoledata.dat"+ext;
          std::ifstream f(filename.c_str());
          if (f) {
            while (!f.eof()) {
              string name;
              int cnt;
              double average;
              double sqsum;
              f >> name >> cnt >> average >> sqsum;
              if(name == "M_x") {
                M_x.reset();
                Average<double> M_xt(average,sqsum,cnt);
                M_x = M_x + M_xt;
              }
              if(name == "M_y") {
                M_y.reset();
                Average<double>M_yt(average,sqsum,cnt);
                M_y = M_y + M_yt;
              }
              if(name == "M_z") {
                M_z.reset();
                Average<double> M_zt(average,sqsum,cnt);
                M_z = M_z + M_zt;
              }
              if(name == "M_x_box") {
                M_x_box.reset();
                Average<double> M_x_boxt(average,sqsum,cnt);
                M_x_box = M_x_box + M_x_boxt;
              }
              if(name == "M_y_box") {
                M_y_box.reset();
                Average<double> M_y_boxt(average,sqsum,cnt);
                M_y_box = M_y_box + M_y_boxt;
              }
              if(name == "M_z_box") {
                M_z_box.reset();
                Average<double> M_z_boxt(average,sqsum,cnt);
                M_z_box = M_z_box + M_z_boxt;
              }
              if(name == "M2") {
                M2.reset();
                Average<double> M2t(average,sqsum,cnt);
                M2 = M2 + M2t;
              }
              if(name == "M2_box") {
                M2_box.reset();
                Average<double> M2_boxt(average,sqsum,cnt);
                M2_box = M2_box + M2_boxt;
              }
              if(name == "diel_std") {
                diel_std.reset();
                Average<double> diel_stdt(average,sqsum,cnt);
                diel_std = diel_std + diel_stdt;
              }
              if(name == "mu_abs") {
                mu_abs.reset();
                Average<double> mu_abst(average,sqsum,cnt);
                mu_abs = mu_abs + mu_abst;
              }
            }
          }
        }

        inline string info() {
          using namespace Faunus::textio;
          std::ostringstream o;
          o << header("Dipole analysis");
          o << indent(SUB) << epsilon_m+subr+"(Tinfoil)" << setw(22) << getDielTinfoil() << ", "+sigma+"=" << diel_std.stdev() << ", "+sigma+"/"+epsilon_m+subr+"=" << (100*diel_std.stdev()/getDielTinfoil()) << percent << endl
          << indent(SUB) << bracket("M"+squared) << setw(27) << pc::eA2D(M2_box.avg(),2) << " Debye"+squared+", "+sigma+"=" << pc::eA2D(M2_box.stdev(),2) << ", "+sigma+"/"+bracket("M"+squared)+"=" << (100*M2_box.stdev()/M2_box.avg()) << percent << endl
          << indent(SUB) << bracket("M") << setw(25) << "( " << pc::eA2D(M_x_box.avg()) << " , " << pc::eA2D(M_y_box.avg()) << " , " << pc::eA2D(M_z_box.avg()) << " ) Debye" << endl 
          << indent(SUBSUB) << sigma << setw(25) << "( " << pc::eA2D(M_x_box.stdev()) << " , " << pc::eA2D(M_y_box.stdev()) << " , " << pc::eA2D(M_z_box.stdev()) << " )" << endl
          << indent(SUB) << bracket("|"+mu+"|") << setw(25) << pc::eA2D(mu_abs.avg()) << " Debye, "+sigma+"=" << pc::eA2D(mu_abs.stdev()) << ", "+sigma+"/"+bracket("|"+mu+"|")+"=" << (100*mu_abs.stdev()/mu_abs.avg()) << percent << endl;
          return o.str();
        }
    };

    /*
     * Perhaps make this a template, taking T=double as parameter?
     */
    template<class T=double>
      class analyzeVector {
        private:
          std::vector<T> noise;
          int N;
          int lag;
          double mu;
          double sigma2;
          double F_X_LB;
          double F_X_BP;
          double F_X_stud;
          double sTd_alpha;
        public:
          inline analyzeVector(std::vector<T> noise_in) {
            noise = noise_in;
            N = noise.size();
            lag = std::round(std::round(std::log(N)));  // More general instead of the maybe more common lag=20'
            lag = N;                                  // Only meaningful for large N (???)
            initMu();
            initVariance();
            F_X_LB = 0;
            F_X_BP = 0;
            F_X_stud = 0;
            sTd_alpha = 0;
          }

          void setLag(int lag_in) {
            lag = lag_in;
          }

          void initMu() {
            mu = 0;
            for(int k = 0;k < N; k++) {
              mu += noise[k];
            }
            mu /= N;
          }

          double getMu() {
            return mu;
          }

          void initVariance() {
            sigma2 = 0;
            for(int k = 0;k < N; k++) {
              sigma2 += noise[k]*noise[k];
            }
            sigma2 -= mu*mu;
          }

          double getVariance() {
            return sigma2;
          }

          // Approximate cumulative distribution function for value alpha in Ljung-Box method
          double getF_X_LB() {
            return F_X_LB;
          }

          // Approximate cumulative distribution function for value alpha in Box-Pierce method
          double getF_X_BP() {
            return F_X_BP;
          }

          // Approximate cumulative distribution function for value alpha in studentTdistribution method
          double getF_X_stud() {
            return F_X_stud;
          }

          // Value of alpha to approximately get F_X_stud in studentTdistribution method (Method fails if (x*x/dof > 1), see line further down)
          double getStdAlpha() {
            return sTd_alpha;
          }

          /**
           * @brief Chi-squared distribution. 
           *        \f$ F(x;k) = \frac{\gamma\left(\frac{k}{2},\frac{x}{2}\right)}{\Gamma\left(\frac{k}{2}\right)} \f$ where \f$ s \f$ is half the lag.
           * 
           * @param significance The significance to reach.
           * @param chi2step Step size in the iteration.
           */
          double getChi2(double significance, double &F_X, double chi2_step=0.01) {
            double x = -chi2_step;
            double level = 1-significance;
            do {
              x += chi2_step;
              F_X = incompleteGamma(x/2)/std::tgamma(double(lag)/2);
            } while (level > F_X);
            return x;
          }

          /**
           * @brief Lower Incomplete Gamma Function.
           *        \f$ \gamma(s,x) = x^se^{-x}\sum_{k=0}^{\infty}\frac{x^k}{s(s+1)...(s+k)} \f$ where \f$ s \f$ is half the lag.
           * 
           * @param x Value to estimate for.
           */
          double incompleteGamma(double x){
            double s = double(lag)/2;
            double sum = 0;
            double term = 1.0/s;
            int n = 1;
            while (term != 0){
              sum += term;
              term *= (x/(s+n++));
            }
            return pow(x,s)*exp(-x)*sum;
          }

          double sampleAutocorrelation(double k) {
            double nom = 0;
            double den = 0;
            for(int t = 0; t < N-k;t++) {
              nom += (noise[t]-mu)*(noise[t+k]-mu);
              den += (noise[t]-mu)*(noise[t]-mu);
            }
            for(int t = N-k; t < N;t++) {
              den += (noise[t]-mu)*(noise[t]-mu);
            }
            return (nom/den);
          }

          /**
           * @brief The Ljung-Box test is a statistical test.
           *
           * This tests if a group of autocorrelations differ from zero, 
           * based on a number of lags. Initially it was developed for ARMA-processes.
           * It is a portmanteau test. More info at DOI: 10.1093/biomet/65.2.297
           *
           * @param alpha The significance of the test. With a probability of \f$ (1-\alpha)*100 \f$ % the result is true.
           */
          bool LjungBox(double alpha) {
            double Q = 0.0;
            F_X_LB = 0;
            for(int k = 0;k < lag; k++)
              Q += sampleAutocorrelation(k)/(N-k);

            if(Q*N*(N+2) > getChi2(alpha,F_X_LB))
              return false;
            return true;
          }

          /**
           * @brief The Box-Pierce test is a statistical test.
           *
           * This tests if a group of autocorrelations differ from zero. 
           * It tests the wider randomness based on a number of lags.
           * This is more simple, and not as accurate, as the Ljung_Box test.
           *
           * @param alpha The significance of the test. With a probability of
           * \f$ (1-\alpha)*100 \f$ % the result is true.
           */
          bool BoxPierce(double alpha) {
            double Q = 0.0;
            F_X_BP = 0;
            for(int k = 0;k < lag; k++)
              Q += sampleAutocorrelation(k);

            if(Q*N > getChi2(alpha,F_X_BP))
              return false;
            return true;
          }

          /**
           * @brief Hypergeometric function.
           *
           * This function uses \f$ x=(x)_1=(x)_2=... \f$ for x=a,x=b and x=c.
           * \f$ F_1(a,b,c;z) = \sum_{n=0}^{\infty}\frac{(a)_n(b)_n}{(c)_n}\frac{z^n}{n!} \f$
           * 
           * @param a Coefficient (usually a vector) 
           * @param b Coefficient (usually a vector) 
           * @param c Coefficient (usually a vector) 
           * @param z Value to estimate for. Only works for \f$ |z|<1 \f$.
           */
          double F2(double a,double b,double c, double z) {
            assert( std::abs(z) < 1 && "|z| is not smaller than 1!");
            int n = 0;
            int nfac = 1;
            double term = 0;
            double sum = 0;
            double coeff = a*b/c;
            do {
              term = coeff*(pow(z,n)/nfac);
              sum += term;
              nfac *= ++n;
            } while(term > 1e-30); // Some cutoff
            return sum;
          }

          /**
           * @brief Student's t-test distribution
           * 
           * @param alpha The significance of the test.
           *        With a probability of \f$ (1-\alpha)*100 \f$ % the result is true.
           * @param dof Degrees of freedom.
           */
          double studentTdistribution(double alpha, int dof, double x_step=0.01) {
            double level = 1 - alpha;
            F_X_stud = 0;
            double x = -x_step;
            do {
              x += 0.01;
              if(x*x/dof > 1) {
                x -= x_step;
                sTd_alpha = 1 - F_X_stud;
                break;
              }
              F_X_stud = 0.5+x*std::tgamma((dof+1)/2)*F2(0.5,(dof+1)/2,1.5,-x*x/dof)/(sqrt(pc::pi*dof)*std::tgamma(dof/2));
            } while (level > F_X_stud);
            return x;
          }

          /**
           * @brief Testing the null hypothesis that the mean is equal to zero. The degrees of freedom used in this test is N − 1. 
           * 
           * @param alpha The significance of the test. With a probability of \f$ (1-\alpha)*100 \f$ % the result is true.
           */
          bool oneSampleT(double alpha) {
            return (studentTdistribution(alpha, N-1) > sqrt(N/sigma2)*mu);
          }

          /**
           * @brief Returns the result of the LjungBox-, BoxPierce- and oneSampleT-tests.
           * 
           * @param alpha The significance of the test. With a probability of \f$ (1-\alpha)*100 \f$ % the result is true.
           */
          std::vector<bool> checkAll(double alpha) {
            std::vector<bool> tests(3);
            tests.at(0) = LjungBox(alpha);
            tests.at(1) = BoxPierce(alpha);
            tests.at(2) = oneSampleT(alpha);
            return tests;
          }

          string info(char w=0) {
            using namespace Faunus::textio;
            std::ostringstream o;
            o << "Sample size: " << N << ", Mean value: " << mu << ", Variance: " << sigma2 << ", Lag: " << lag << endl;
            o << "F_X(Ljung-Box): " << F_X_LB << ", F_X(Box_Pierce): " << F_X_BP << ", Significance(oneSampleT): " << sTd_alpha << endl;
            return o.str();
          }
      };
  }//namespace
}//namespace
#endif
