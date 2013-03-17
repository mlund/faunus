#ifndef faunus_energy_h
#define faunus_energy_h

#ifndef SWIG
#include <faunus/common.h>
#include <faunus/point.h>
#include <faunus/group.h>
#include <faunus/inputfile.h>
#include <faunus/space.h>
#include <faunus/textio.h>
#include <faunus/potentials.h>
#include <faunus/auxiliary.h>
#include <faunus/analysis.h>
#endif

// http://publib.boulder.ibm.com/infocenter/iadthelp/v8r0/index.jsp?topic=/com.ibm.xlcpp111.linux.doc/language_ref/variadic_templates.html

namespace Faunus {

  /**
   * @brief Classes/templates that calculate the energy of particles, groups, system.
   *
   * This namespace containes classes and templates for calculating energies of the
   * system. This can be pair interactions, external potentials, long range corrections,
   * constraints etc.
   * All energy calculation classes inherit from the base class Energy::Energybase and
   * thus have a common interface. Of special significance is the class
   * Energy::Hamiltonian that can sum the energies from an arbitrary number of base
   * classes and thus act as a true Hamiltonian.
   */
  namespace Energy {


    /**
     *  @brief Base class for energy evaluation
     *
     *  This base class defines functions for evaluating interactions between particles,
     *  groups, external potentials etc. By default all energy functions return ZERO
     *  and derived classes are expected only to implement functions relevant for certain
     *  properties. I.e. a derived class for non-bonded interactions are not expected to
     *  implement `i_internal()`, for example.
     *
     *  @note All energy functions are expected to return energies in units of kT.
     *  @todo Add setVolume() function such that each derived class may have it's own
     *        Geometry instance (if needed). This will significantly simplify the
     *        Hamiltonian class and increase performance by avoiding calling
     *        distance functions via the "geo" pointer.
     */
    class Energybase {
      private:
        virtual string _info()=0;
        char w; //!< Width of info output
      protected:
        Geometry::Geometrybase* geo; //!< Pointer to geometry used to calculate interactions
      public:
        string name;                                          //!< Short informative name
        Energybase();                                         //!< Constructor
        virtual ~Energybase();                                //!< Destructor
        virtual Geometry::Geometrybase& getGeometry();        //!< Geometry used for interactions
        bool setGeometry( Geometry::Geometrybase& );          //!< Set Geometrybase
        virtual void setTemperature(double);                  //!< Set temperature for interactions
        virtual void setVolume(double);                       //!< Set volume of used Geometry
        virtual double p2p(const particle&, const particle&); // Particle-particle energy
        virtual Point f_p2p(const particle&, const particle&); // Particle-particle force
        virtual double all2p(const p_vec&, const particle&);  // Particle-Particle vector energy
        virtual double all2all(const p_vec&);                 // All inter-particle energies (N^2)
        virtual double i2i(const p_vec&, int, int);           // i'th particle with j'th particle
        virtual double i2g(const p_vec&, Group &, int);       // i'th particle with group
        virtual double i2all(const p_vec&, int);              // i'th particle with all other particles
        virtual double i_external(const p_vec&, int);         // internal energy of i'th particle
        virtual double i_internal(const p_vec&, int);         // External energy of i'th particle
        virtual double p_external(const particle&);           // External energy of particle
        double i_total(const p_vec&, int);                    // Total energy of i'th particle = i2all + i_external + i_internal
        virtual double g2g(const p_vec&, Group&, Group&);     // Group-Group energy
        virtual double g2all(const p_vec&, Group&);           // Energy of Group with all other particles
        virtual double g_external(const p_vec&, Group&);      // External energy of group
        virtual double g_internal(const p_vec&, Group&);      // Internal energy of group
        virtual double v2v(const p_vec&, const p_vec&);       // Particle vector-Particle vector energy
        virtual double external();                            // External energy - pressure, for example.
        virtual string info();                                //!< Information

        virtual void field(const p_vec&, std::vector<Point>&);//!< Calculate electric field on all particles

        virtual void trialUpdate(const Space&, std::set<int>&) {};
        virtual void acceptUpdate() {};
        virtual void rejectUpdate() {};
    };

    /**
     * @brief Energy class for non-bonded interactions.
     *
     * `Tpairpot` is expected to be a pair potential with the following properties:
     *
     * - `Tpairpot(const InputMap&)`
     * - `double Tpairpot::operator()(const particle&, const particle&, double sqdist))`
     *
     * For a list of implemented potentials, see the Faunus::Potential namespace.
     */
    template<class Tpairpot, class Tgeometry>
      class Nonbonded : public Energybase {
        protected:
          string _info() {
            return pairpot.info(25);
          }
        public:
          Tgeometry geometry;
          Tpairpot pairpot;
          Nonbonded(InputMap &in) : geometry(in), pairpot(in) {
            static_assert(
                std::is_base_of<Potential::PairPotentialBase,Tpairpot>::value,
                "Tpairpot must be a pair potential" );
            name="Nonbonded N" + textio::squared + " - " + pairpot.name;
            setGeometry(geometry);
          }

          Geometry::Geometrybase& getGeometry() {
            geo=&geometry;
            return Energybase::getGeometry();
          }

          //!< Particle-particle energy (kT)
          inline double p2p(const particle &a, const particle &b) FOVERRIDE {
            return pairpot( a,b,geometry.sqdist(a,b));
          }

          //!< Particle-particle force (kT/Angstrom)
          inline Point f_p2p(const particle &a, const particle &b) FOVERRIDE {
            return pairpot.force( a,b,geometry.sqdist(a,b),geometry.vdist(a,b));
          }

          double all2p(const p_vec &p, const particle &a) FOVERRIDE {
            double u=0;
            for (auto &b : p)
              u+=pairpot(a,b,geometry.sqdist(a,b));
            return u;
          }

          double all2all(const p_vec &p) FOVERRIDE {
            int n=p.size();
            double u=0;
            for (int i=0; i<n-1; ++i)
              for (int j=i+1; j<n; ++j)
                u+=pairpot( p[i],p[j],geometry.sqdist(p[i],p[j]) );
            return u;
          }

          double i2i(const p_vec &p, int i, int j) FOVERRIDE {
            return pairpot( p[i], p[j], geometry.sqdist( p[i], p[j]) );
          }

          double i2g(const p_vec &p, Group &g, int j) FOVERRIDE {
            double u=0;
            if ( !g.empty() ) {
              int len=g.back()+1;
              if ( g.find(j) ) {   //j is inside g - avoid self interaction
                for (int i=g.front(); i<j; i++)
                  u+=pairpot(p[i],p[j],geometry.sqdist(p[i],p[j]));
                for (int i=j+1; i<len; i++)
                  u+=pairpot(p[i],p[j],geometry.sqdist(p[i],p[j]));
              } else              //simple - j not in g
                for (int i=g.front(); i<len; i++)
                  u+=pairpot( p[i], p[j], geometry.sqdist(p[i],p[j]));
            }
            return u;  
          }

          double i2all(const p_vec &p, int i) FOVERRIDE {
            assert(i>=0 && i<int(p.size()) && "index i outside particle vector");
            double u=0;
            int n=(int)p.size();
            for (int j=0; j<i; ++j)
              u+=pairpot( p[i], p[j], geometry.sqdist(p[i],p[j]) );
            for (int j=i+1; j<n; ++j)
              u+=pairpot( p[i], p[j], geometry.sqdist(p[i],p[j]) );
            return u;
          }

          double g2g(const p_vec &p, Group &g1, Group &g2) FOVERRIDE {
            double u=0;
            if (!g1.empty())
              if (!g2.empty()) {
                // IN CASE ONE GROUP IS A SUBGROUP OF THE OTHER
                if (g1.find(g2.front()))
                  if (g1.find(g2.back())) {  // g2 is a subgroup of g1
                    assert(g1.size()>=g2.size());
                    for (int i=g1.front(); i<g2.front(); i++)
                      for (auto j : g2)
                        u+=pairpot(p[i],p[j],geometry.sqdist(p[i],p[j]));
                    for (int i=g2.back()+1; i<=g1.back(); i++)
                      for (auto j : g2)
                        u+=pairpot(p[i],p[j],geometry.sqdist(p[i],p[j]));
                    return u;
                  }
                if (g2.find(g1.front()))
                  if (g2.find(g1.back())) {  // g1 is a subgroup of g2
                    assert(g2.size()>=g1.size());
                    for (int i=g2.front(); i<g1.front(); i++)
                      for (auto j : g1)
                        u+=pairpot(p[i],p[j],geometry.sqdist(p[i],p[j]));
                    for (int i=g1.back()+1; i<=g2.back(); i++)
                      for (auto j : g1)
                        u+=pairpot(p[i],p[j],geometry.sqdist(p[i],p[j]));
                    return u;
                  }

                // IN CASE BOTH GROUPS ARE INDEPENDENT (DEFAULT)
                int ilen=g1.back()+1, jlen=g2.back()+1;
#pragma omp parallel for reduction (+:u) schedule (dynamic)
                for (int i=g1.front(); i<ilen; ++i)
                  for (int j=g2.front(); j<jlen; ++j)
                    u+=pairpot(p[i],p[j],geometry.sqdist(p[i],p[j]));
              }
            return u;
          }

          double g2all(const p_vec &p, Group &g) FOVERRIDE {
            double u=0;
            if (!g.empty()) {
              int ng=g.back()+1, np=p.size();
#pragma omp parallel for reduction (+:u)
              for (int i=g.front(); i<ng; ++i) {
                for (int j=0; j<g.front(); j++)
                  u += pairpot( p[i], p[j], geometry.sqdist(p[i],p[j]) );
                for (int j=ng; j<np; j++)
                  u += pairpot( p[i], p[j], geometry.sqdist(p[i],p[j]) );
              }
            }
            return u;
          }

          double g_internal(const p_vec &p, Group &g) FOVERRIDE { 
            double u=0;
            const int b=g.back(), f=g.front();
            if (!g.empty())
              for (int i=f; i<b; ++i)
                for (int j=i+1; j<=b; ++j)
                  u+=pairpot(p[i],p[j],geometry.sqdist(p[i],p[j]));
            return u;
          }

          double v2v(const p_vec &p1, const p_vec &p2) FOVERRIDE {
            double u=0;
            for (auto &i : p1)
              for (auto &j : p2)
                u+=p2p(i,j);
            return u;
          }
      };

    /**
     * @brief Energy class for non-bonded interactions.
     *
     * `Tpairpot` is expected to be a pair potential with the following properties:
     *
     * - `Tpairpot(const InputMap&)`
     * - `double Tpairpot::operator()(const particle&, const particle&, double sqdist))`
     *
     * For a list of implemented potentials, see the Faunus::Potential namespace.
     */
    template<class Tpairpot, class Tgeometry>
      class NonbondedVector : public Energybase {
        protected:
          string _info() {
            return pairpot.info(25);
          }
        public:
          Tgeometry geometry;
          Tpairpot pairpot;
          NonbondedVector(InputMap &in) : geometry(in), pairpot(in) {
            name="Nonbonded N" + textio::squared + " - " + pairpot.name;
            geo=&geometry;
          }

          Geometry::Geometrybase& getGeometry() {
            geo=&geometry;
            return Energybase::getGeometry();
          }

          //!< Particle-particle energy (kT)
          inline double p2p(const particle &a, const particle &b) FOVERRIDE {
            return pairpot( a,b,geometry.vdist(a,b));
          }

          //!< Particle-particle force (kT/Å)
          inline Point f_p2p(const particle &a, const particle &b) FOVERRIDE {
            return pairpot.force( a,b,geometry.sqdist(a,b),geometry.vdist(a,b));
          }

          double all2p(const p_vec &p, const particle &a) FOVERRIDE {
            double u=0;
            for (auto &b : p)
              u+=pairpot(a,b,geometry.vdist(a,b));
            return u;
          }

          double all2all(const p_vec &p) FOVERRIDE {
            int n=p.size();
            double u=0;
            for (int i=0; i<n-1; ++i)
              for (int j=i+1; j<n; ++j)
                u+=pairpot( p[i],p[j],geometry.vdist(p[i],p[j]) );
            return u;
          }

          double i2i(const p_vec &p, int i, int j) FOVERRIDE {
            return pairpot( p[i], p[j], geometry.vdist( p[i], p[j]) );
          }

          double i2g(const p_vec &p, Group &g, int j) FOVERRIDE {
            double u=0;
            if ( !g.empty() ) {
              int len=g.back()+1;
              if ( g.find(j) ) {   //j is inside g - avoid self interaction
                for (int i=g.front(); i<j; i++)
                  u+=pairpot(p[i],p[j],geometry.vdist(p[i],p[j]));
                for (int i=j+1; i<len; i++)
                  u+=pairpot(p[i],p[j],geometry.vdist(p[i],p[j]));
              } else              //simple - j not in g
                for (int i=g.front(); i<len; i++)
                  u+=pairpot( p[i], p[j], geometry.vdist(p[i],p[j]));
            }
            return u;
          }

          double i2all(const p_vec &p, int i) FOVERRIDE {
            assert(i>=0 && i<int(p.size()) && "index i outside particle vector");
            double u=0;
            int n=(int)p.size();
            for (int j=0; j<i; ++j)
              u+=pairpot( p[i], p[j], geometry.vdist(p[i],p[j]) );
            for (int j=i+1; j<n; ++j)
              u+=pairpot( p[i], p[j], geometry.vdist(p[i],p[j]) );
            return u;
          }

          double g2g(const p_vec &p, Group &g1, Group &g2) FOVERRIDE {
            double u=0;
            if (!g1.empty())
              if (!g2.empty()) {
                // IN CASE ONE GROUP IS A SUBGROUP OF THE OTHER
                if (g1.find(g2.front()))
                  if (g1.find(g2.back())) {  // g2 is a subgroup of g1
                    assert(g1.size()>=g2.size());
                    for (int i=g1.front(); i<g2.front(); i++)
                      for (auto j : g2)
                        u+=pairpot(p[i],p[j],geometry.vdist(p[i],p[j]));
                    for (int i=g2.back()+1; i<=g1.back(); i++)
                      for (auto j : g2)
                        u+=pairpot(p[i],p[j],geometry.vdist(p[i],p[j]));
                    return u;
                  }
                if (g2.find(g1.front()))
                  if (g2.find(g1.back())) {  // g1 is a subgroup of g2
                    assert(g2.size()>=g1.size());
                    for (int i=g2.front(); i<g1.front(); i++)
                      for (auto j : g1)
                        u+=pairpot(p[i],p[j],geometry.vdist(p[i],p[j]));
                    for (int i=g1.back()+1; i<=g2.back(); i++)
                      for (auto j : g1)
                        u+=pairpot(p[i],p[j],geometry.vdist(p[i],p[j]));
                    return u;
                  }

                // IN CASE BOTH GROUPS ARE INDEPENDENT (DEFAULT)
                int ilen=g1.back()+1, jlen=g2.back()+1;
#pragma omp parallel for reduction (+:u) schedule (dynamic)
                for (int i=g1.front(); i<ilen; ++i)
                  for (int j=g2.front(); j<jlen; ++j)
                    u+=pairpot(p[i],p[j],geometry.vdist(p[i],p[j]));
              }
            return u;
          }

          double g2all(const p_vec &p, Group &g) FOVERRIDE {
            double u=0;
            if (!g.empty()) {
              int ng=g.back()+1, np=p.size();
#pragma omp parallel for reduction (+:u)
              for (int i=g.front(); i<ng; ++i) {
                for (int j=0; j<g.front(); j++)
                  u += pairpot( p[i], p[j], geometry.vdist(p[i],p[j]) );
                for (int j=ng; j<np; j++)
                  u += pairpot( p[i], p[j], geometry.vdist(p[i],p[j]) );
              }
            }
            return u;
          }

          double g_internal(const p_vec &p, Group &g) FOVERRIDE {
            double u=0;
            if (!g.empty()) {
              int step=1,n=g.back()+1;
              for (int i=g.front(); i<n-step; i++)
                for (int j=g.front()+step*((i-g.front())/step+1); j<n; j++)
                  u+=pairpot(p[i],p[j],geometry.vdist(p[i],p[j]));
            }
            return u;
          }

          double v2v(const p_vec &p1, const p_vec &p2) FOVERRIDE {
            double u=0;
            for (auto &i : p1)
              for (auto &j : p2)
                u+=p2p(i,j);
            return u;
          }

          /**
           * Calculates the electric field on all particles
           * and stores (add) in the vector `E`.
           *
           * @param p Particle vector
           * @param E Holds field on each particle. Must be of size N.
           */
          void field(const p_vec &p, std::vector<Point> &E) {
            assert(p.size()==E.size());
            size_t i=0;
            for (auto &pi : p) {
              for (auto &pj : p)
                if (&pi!=&pj)
                  E[i]+=pairpot.field(pj, geometry.vdist(pi,pj));
              i++;
            }
          }
      };

    /**
     * @brief Nonbonded interactions with group-group cutoff
     *
     * This class re-implements the `g2g` energy function (group to group)
     * and checks of the mass center distance is larger than a certain
     * cut-off. If so zero is returned - otherwise the usual `g2g` function
     * from the parent class. The cut-off is observed only between groups
     * where `isMolecular()` is true.
     */
    template<class Tpairpot, class Tgeometry>
      class NonbondedCut : public Nonbonded<Tpairpot, Tgeometry> {
        using Nonbonded<Tpairpot,Tgeometry>::geo;
        using Nonbonded<Tpairpot,Tgeometry>::name;
        private:
        Space* spcPtr;
        unsigned int cnt, cntfull;
        double sqcut; //!< Squared cut-off mass-center distance
        string _info() {
          char w=25;
          using namespace textio;
          std::ostringstream o;
          o << Nonbonded<Tpairpot, Tgeometry>::_info()
            << pad(SUB,w,"Group-Group cut-off") << sqrt(sqcut) << _angstrom << "\n";
          if (cnt>0)
            o << pad(SUB,w,"Cut-off fraction") << 1-cntfull/double(cnt) << "\n";
          return o.str();
        }
        public:
        NonbondedCut(InputMap &in) : Nonbonded<Tpairpot, Tgeometry>(in) {
          sqcut = pow( in.get<double>("g2g_cutoff",pc::infty), 2);
          name+="(Molecular Group Cutoff)";
          spcPtr=nullptr;
          cnt=cntfull=0;
        }

        /**
         * @brief Specify simulation Space. Needed to distinguish trial configurations
         *        from old ones.
         */
        void setSpace(Space &spc) { spcPtr=&spc; }

        double g2g(const p_vec &p, Group &g1, Group &g2) {
          cnt++;
          assert(spcPtr!=nullptr && "You forgot to set Space.");
          if (g1.isMolecular()) {         // only between molecular groups
            if (g2.isMolecular()) {       // -//-
              if ( &p==&spcPtr->trial ) { // new or old particle vector?
                if ( geo->sqdist(g1.cm_trial, g2.cm_trial) > sqcut )
                  return 0;
              } else
                if ( geo->sqdist(g1.cm, g2.cm) > sqcut )
                  return 0;
            }
          }
          cntfull++;
          return Nonbonded<Tpairpot, Tgeometry>::g2g(p,g1,g2); // full energy
        }
      };

    /**
     * @brief Energy class for hard-sphere overlap.
     */
    template<class Tgeometry>
      class HardSphereOverlap : public Energybase {
        private:
          Tgeometry geo;
          Potential::HardSphere hs;
        public:
          HardSphereOverlap(InputMap &mcp) : geo(mcp) {} 
          inline double i2i(const p_vec &p, int i, int j) FOVERRIDE {
            return hs(p[i], p[j], geo.sqdist(p[i], p[j]) );
          }
          double all2all(const p_vec &p) FOVERRIDE {
            for (auto i=p.begin(); i!=p.end()-1; ++i)
              for (auto j=i+1; j!=p.end(); ++j)
                if ( hs(*i,*j, geo.sqdist(*i,&j) )==pc::infty )
                  return pc::infty;
            return 0;
          }
          double g2all(const p_vec &p, Group &g) FOVERRIDE {
            if (g.empty())
              return 0;
            for (int i=g.front(); i<=g.back(); i++) {
              for (int j=0; j<g.front(); j++)
                if ( i2i(p,i,j)==pc::infty )
                  return pc::infty;
              for (int j=g.back()+1; j<(int)p.size(); j++)
                if ( i2i(p,i,j)==pc::infty )
                  return pc::infty;
              return 0;
            }
          }
          double g2g(const p_vec &p, Group &g1, Group &g2) FOVERRIDE {
            for (int i=g1.front(); i<=g1.back(); ++i)
              for (int j=g2.front(); j<=g2.back(); ++j)
                if ( i2i(p,i,j)==pc::infty )
                  return pc::infty;
            return 0;
          }
      };

    /**
     * @brief Class for handling bond pairs
     *
     * Takes care of bonded interactions and can handle mixed bond types. If you create bond BETWEEN
     * groups, make sure to set the `CrossGroupBonds` to `true`.
     *
     * Example:
     *
     *     vector<particle> p(...);            // particle vector
     *     int i=10, j=11;                     // particle index
     *     Energy::Bonded b;
     *     b.add(i, j, Potential::Harmonic(0.1,5.0) );
     *     std::cout << b.info();
     *     double rij2 = ... ;                 // squared distance between i and j
     *     double u = b(i,j)( p[i], p[j], rij2 ); // i j bond energy in kT
     *
     * @date Lund, 2011-2012
     */
    class Bonded : public Energybase, public pair_list<Potential::PairPotentialBase> {
      private:
        string _info();
      public:
        Bonded();
        Bonded(Geometry::Geometrybase&);
        double i2i(const p_vec&, int, int) FOVERRIDE;      //!< Bond energy i with j
        double i2all(const p_vec&, int) FOVERRIDE;         //!< All bonds w. i'th particle
        double g_internal(const p_vec&, Group&) FOVERRIDE; //!< Internal bonds in Group, only
        double g2g(const p_vec&, Group&, Group&) FOVERRIDE;//!< Bonds between groups
        double total(const p_vec&);                        //!< Sum all known bond energies
        bool CrossGroupBonds; //!< Set to true if bonds cross groups (slower!). Default: false
    };

    /**
     * @brief Energy from external pressure for use in the NPT-ensemble.
     *
     * @details This will count the number of particles in the system and
     * calculate the energy contribution to the pressure at a given
     * volume and external pressure.
     *
     * The system energy is
     * \f$\beta u = \beta pV - \ln V - N\ln V\f$
     * where the first two terms are returned by `external()` while the last
     * term is obtained by summing `g_external()` over molecular
     * and or atomic groups.
     * If applied on an atomic group, `N` will be set to the number of
     * atoms in the group, while for a molecular group `N=1`.
     *
     * @date Lund, 2011
     */
    class ExternalPressure : public Energy::Energybase {
      private:
        double P; //!< Pressure, p/kT
        string _info();
      public:
        ExternalPressure(Geometry::Geometrybase&, double);
        double external() FOVERRIDE;  //!< External energy working on system. pV/kT-lnV
        double g_external(const p_vec&, Group&) FOVERRIDE; //!< External energy working on group
    };

    /**
     * @brief External energy that will keep specific groups in a sub-volume of the system
     *
     * This energy class will check if particles in specific groups are located within a
     * rectangular box, spanned by two vector points, `upper` and `lower`. If outside
     * an infinite energy is returned. This is useful for constraining molecules in specific
     * parts of the simulation container.
     * Derived classes can re-implement the virtual
     * `outside()` function which should return `true` if a given point falls outside the
     * allowed region. Note that the current implementation can be problematic with
     * containers with periodic boundaries as `outside()` uses absolute positions.
     *
     * Example:
     *
     *     Energy::Hamiltonian pot;
     *     ...
     *     auto restricted = pot.create( Energy::RestrictedVolume(imap) );
     *     restricted->groups.push_back( &mygroup );
     *
     * @date Lund, 2012
     */
    class RestrictedVolume : public Energy::Energybase {
      private:
        Point upper, lower;
        string _info();
      protected:
        virtual bool outside(const Point&); //!< Determines if particle is outside allowed region
      public:
        std::vector<Group*> groups;              //!< List of groups to restrict
        RestrictedVolume(InputMap&, string="vconstrain"); //!< Constructor
        double g_external(const p_vec&, Group&) FOVERRIDE; //!< External energy working on group
    };

    /**
     * @brief As Energy::RestrictedVolume but restrictions are applied only on the mass center
     * instead of all particles in group.
     */
    class RestrictedVolumeCM : public Energy::RestrictedVolume {
      public:
        RestrictedVolumeCM(InputMap&, string="vconstrain"); //!< Constructor
        double g_external(const p_vec&, Group&) FOVERRIDE; //!< External energy working on group
        double i_external(const p_vec&, int);              //!< External energy working on single particle
    };

    /**
     * @brief Collection of Energybases that when summed give the Hamiltonian
     *
     * This class is used to collect several Energybase derivatives into a full hamiltonian.
     * The following example demonstrated how one can generate a Hamiltonian for bonded as
     * well as non-bonded interactions:
     * \code
     * Energy::Hamiltonian pot;
     * auto nbPtr = pot.create( Energy::Nonbonded<Tpairpot,Tgeometry>(in) );
     * auto bPtr = pot.create( Energy::Bonded() );
     * cout << pot.info();
     * \endcode
     *
     * Notice that we do not need to specify a Geometry for the Bonded energy class as this information
     * is simply passed on from the first added potential.
     *
     * \author Mikael Lund
     * @date Lund, 2011
     */
    class Hamiltonian : public Energybase {
      typedef shared_ptr<Energybase> baseptr;
      private:
      vector<baseptr> created;      //!< smart pointer list of *created* energy classes
      string _info();
      vector<Energybase*> baselist; //!< Pointer list to energy classes to be summed
      public:
      Hamiltonian();
      void setVolume(double);       //!< Set volume of all contained energy classes
      void setTemperature(double);  //!< Set temperature of all contained energy classes

      /**
       * @brief Create and add an energy class to energy list
       */
      template<typename Tenergychild> shared_ptr<Tenergychild> create(Tenergychild c) {
        shared_ptr<Tenergychild> childptr( new Tenergychild(c) );
        childptr->getGeometry(); // not pretty...need to update geo pointer for i.e. nonbonded class
        created.push_back(childptr);
        add(*childptr);
        return childptr;
      }

      void add(Energybase&); //!< Add existing energy class to list
      double p2p(const particle&, const particle&) FOVERRIDE;
      Point f_p2p(const particle&, const particle&) FOVERRIDE;
      double all2p(const p_vec&, const particle&) FOVERRIDE;
      double all2all(const p_vec&) FOVERRIDE;
      double i2i(const p_vec&, int, int) FOVERRIDE;
      double i2g(const p_vec&, Group&, int) FOVERRIDE;
      double i2all(const p_vec&, int) FOVERRIDE;
      double i_external(const p_vec&, int) FOVERRIDE;
      double i_internal(const p_vec&, int) FOVERRIDE;
      double g2g(const p_vec&, Group&, Group&) FOVERRIDE;
      double g2all(const p_vec&, Group&) FOVERRIDE;
      double g_external(const p_vec&, Group&) FOVERRIDE;
      double g_internal(const p_vec&, Group&) FOVERRIDE;
      double external() FOVERRIDE;
      double v2v(const p_vec&, const p_vec&) FOVERRIDE;
      void field(const p_vec&, std::vector<Point>&) FOVERRIDE;
    };

    /**
     * @brief Constrain two group mass centra within a certain distance interval [mindist:maxdist]
     * \author Mikael Lund
     * @date Lund, 2012
     * \todo Prettify output
     *
     * This energy class will constrain the mass center separation between selected groups to a certain
     * interval. This can be useful to sample rare events and the constraint is implemented as an external
     * group energy that return infinity if the mass center separation are outside the defined range.
     * An arbitrary number of group pairs can be added with the addPair() command, although one would
     * rarely want to have more than one.
     * In the following example,
     * the distance between \c mygroup1 and \c mygroup2 are constrained to the range \c [10:50] angstrom:
     * \code
     * Energy::Hamiltonian pot;
     * auto nonbonded = pot.create( Energy::Nonbonded<Tpairpot,Tgeometry>(mcp) );
     * auto constrain = pot.create( Energy::MassCenterConstrain(pot.getGeometry()) );
     * constrain->addPair( mygroup1, mygroup2, 10, 50); 
     * \endcode
     */
    class MassCenterConstrain : public Energy::Energybase {
      private:
        string _info();
        struct data {
          double mindist, maxdist;
        };
        std::map< opair<Faunus::Group*>, data> gmap;
      public:
        MassCenterConstrain(Geometry::Geometrybase&);      //!< Constructor
        void addPair(Group&, Group&, double, double);      //!< Add constraint between two groups
        double g_external(const p_vec&, Group&) FOVERRIDE; //!< Constrain treated as external potential
    };

    /**
     * @brief Dummy energy class that sums missed energy changes to avoid energy drifts
     *
     * This energy function is designed to be used with Move::Movebase classes
     * that returns energy changes not detectable in the energy drift checkup
     * routines. The idea is simply to sum the energy change discrepancy and treat
     * this as an external potential. Use together with Energy::Hamiltonian.
     */
    class EnergyRest : public Energy::Energybase {
      private:
        double usum;
        string _info();
      public:
        EnergyRest();
        void add(double du); //!< Add energy change disrepancy, dU = U(metropolis) - U(as in drift calculation)
        double external() FOVERRIDE;  //!< Dumme rest treated as external potential to whole system
    };

    /**
     * @brief Sum additive external potential on particles.
     * @tparam Texpot External potential typically derived from
     *         `Potential::ExternalPotentialBase`
     */
    template<typename Texpot>
      class ExternalPotential : public Energybase {
        private:
          string _info() {
            return expot.info();
          }
        public:
          Texpot expot;
          ExternalPotential(InputMap &in) : expot(in) {
            name="External Potential ("+expot.name+")";
          }
          double p_external(const particle &p) FOVERRIDE {
            return expot(p); 
          }
          double i_external(const p_vec &p, int i) FOVERRIDE {
            return p_external( p[i] );
          }
          double g_external(const p_vec &p, Group &g) FOVERRIDE {
            double u=0;
            for (auto i : g)
              u+=p_external(p[i]);
            return u;
          }
      };

    /**
     * @brief Calculates the total system energy
     *
     * For a given particle vector, space, and energy class we try to calculate the
     * total energy taking into account inter- and intra-molecular interactions as well
     * as external potentials. While this may not work for all systems it may be a useful
     * first guess. This is the default energy routine for Move::ParallelTempering and may
     * also be used for checking energy drifts.
     */
    double systemEnergy(Space&, Energy::Energybase&, const p_vec&);

    /* typedefs */

    // template aliasing requires gcc 4.7+
    // template<class Tpairpot, class Tgeometry>
    //   using PairWise = Nonbonded<Tpairpot, Tgeometry>;

  }//Energy namespace
}//Faunus namespace
#endif
