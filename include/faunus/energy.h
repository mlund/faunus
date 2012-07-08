#ifndef faunus_energy_h
#define faunus_energy_h

#include <faunus/common.h>
#include <faunus/point.h>
#include <faunus/group.h>
#include <faunus/textio.h>
#include <faunus/potentials.h>
#include <faunus/auxiliary.h>

// http://publib.boulder.ibm.com/infocenter/iadthelp/v8r0/index.jsp?topic=/com.ibm.xlcpp111.linux.doc/language_ref/variadic_templates.html
//
//
using namespace std;

namespace Faunus {

  /*!
   * \brief Classes related to energy evaluation, pair potentials, bonds etc.
   */
  namespace Energy {


    /*!
     *  \brief Base class for energy evaluation
     *
     *  This base class defines functions for evaluating interactions between particles,
     *  groups, external potentials etc. By default all energy functions return ZERO
     *  and derived classes are expected only to implement functions relevant for certain
     *  properties. I.e. a derived class for non-bonded interactions are not expected to
     *  implement i_internal(), for example.
     *
     *  \note All energy functions are expected to return energies in units of kT.
     *  \todo Add setVolume() function such that each derived class may have it's own
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
        string name;                                          // Short informative name
        Energybase();
        virtual ~Energybase();
        virtual Geometry::Geometrybase& getGeometry();        // Reference to geometry used for interactions
        bool setGeometry( Geometry::Geometrybase& );          // Set Geometrybase
        virtual void setTemperature(double);                  //!< Set temperature for interactions
        virtual void setVolume(double);                       //!< Set volume of used Geometry
        virtual double p2p(const particle&, const particle&); // Particle-particle energy
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
    };

    /*!
     * \brief Energy class for non-bonded interactions.
     *
     * Tpotential is expected to be a pair potential with the following
     * properties:
     * \li pair(const InputMap&);
     * \li double pairpotconst particle&, const particle&, double sqdist);
     * \li double pait.tokT();
     */
    template<class Tpairpot, class Tgeometry>
      class Nonbonded : public Energybase {
        private:
          string _info() {
            return pairpot.info(25);
          }
        public:
          Tgeometry geometry;
          Tpairpot pairpot;
          Nonbonded(InputMap &in) : geometry(in), pairpot(in) {
            name="Nonbonded N" + textio::squared + " - " + pairpot.name;
            geo=&geometry;
          }

          Geometry::Geometrybase& getGeometry() {
            geo=&geometry;
            return Energybase::getGeometry();
          }

          //!< Particle-particle energy (kT)
          inline double p2p(const particle &a, const particle &b) FOVERRIDE {
            return pairpot( a,b,geometry.sqdist(a,b))*pairpot.tokT();
          }

          double all2p(const p_vec &p, const particle &a) FOVERRIDE {
            double u=0;
            for (auto &b : p)
              u+=pairpot(a,b,geometry.sqdist(a,b));
            return u*pairpot.tokT();
          }

          double all2all(const p_vec &p) FOVERRIDE {
            int n=p.size();
            double u=0;
            for (int i=0; i<n-1; ++i)
              for (int j=i+1; j<n; ++j)
                u+=pairpot( p[i],p[j],geometry.sqdist(p[i],p[j]) );
            return u*pairpot.tokT();
          }

          double i2i(const p_vec &p, int i, int j) FOVERRIDE {
            return pairpot.tokT() * pairpot( p[i], p[j], geometry.sqdist( p[i], p[j]) );
          }

          double i2g(const p_vec &p, Group &g, int j) FOVERRIDE {
            if (g.empty())
              return 0;
            double u=0;
            int len=g.back()+1;
            if (j>=g.front() && j<=g.back()) {   //j is inside g - avoid self interaction
              for (int i=g.front(); i<j; i++)
                u+=pairpot(p[i],p[j],geometry.sqdist(p[i],p[j]));
              for (int i=j+1; i<len; i++)
                u+=pairpot(p[i],p[j],geometry.sqdist(p[i],p[j]));
            } else                        //simple - j not in g
              for (int i=g.front(); i<len; i++)
                u+=pairpot( p[i], p[j], geometry.sqdist(p[i],p[j]));
            return pairpot.tokT()*u;  
          }

          double i2all(const p_vec &p, int i) FOVERRIDE {
            double u=0;
            int n=(int)p.size();
            for (int j=0; j<i; ++j)
              u+=pairpot( p[i], p[j], geometry.sqdist(p[i],p[j]) );
            for (int j=i+1; j<n; ++j)
              u+=pairpot( p[i], p[j], geometry.sqdist(p[i],p[j]) );
            return u*pairpot.tokT();
          }

          double g2g(const p_vec &p, Group &g1, Group &g2) FOVERRIDE {
            double u=0;
            if (g1.empty() || g2.empty())
              return u;
            int ilen=g1.back()+1, jlen=g2.back()+1;
#pragma omp parallel for reduction (+:u) schedule (dynamic)
            for (int i=g1.front(); i<ilen; ++i)
              for (int j=g2.front(); j<jlen; ++j)
                u+=pairpot(p[i],p[j],geometry.sqdist(p[i],p[j]));
            return pairpot.tokT()*u;
          }

          double g2all(const p_vec &p, Group &g) FOVERRIDE {
            double u=0;
            if (g.empty())
              return u;
            int ng=g.back()+1, np=p.size();
#pragma omp parallel for reduction (+:u)
            for (int i=g.front(); i<ng; ++i) {
              for (int j=0; j<g.front(); j++)
                u += pairpot( p[i], p[j], geometry.sqdist(p[i],p[j]) );
              for (int j=ng; j<np; j++)
                u += pairpot( p[i], p[j], geometry.sqdist(p[i],p[j]) );
            }
            return u*pairpot.tokT();
          }

          double g_internal(const p_vec &p, Group &g) FOVERRIDE { 
            if (g.empty())
              return 0;
            double u=0;
            int step=1,n=g.back()+1;
            for (int i=g.front(); i<n-step; i++)
              for (int j=g.front()+step*((i-g.front())/step+1); j<n; j++)
                u+=pairpot(p[i],p[j],geometry.sqdist(p[i],p[j]));
            return pairpot.tokT()*u;
          }

          double v2v(const p_vec &p1, const p_vec &p2) FOVERRIDE {
            double u=0;
            for (auto &i : p1)
              for (auto &j : p2)
                u+=p2p(i,j);
            return u;
          }
      };

    template<class Tpairpot, class Tgeometry>
      class Nonbonded_CG : public Nonbonded<Tpairpot, Tgeometry> {
        using Nonbonded<Tpairpot,Tgeometry>::geo;
        using Nonbonded<Tpairpot,Tgeometry>::name;
        public:
        double cut;
        Nonbonded_CG(InputMap &in) : Nonbonded<Tpairpot, Tgeometry>(in) {
          name+="(Molecular Group CG)";
        }

        /*!
         * \warning Why sqdist to trial?? Are we sure this g2g is called? (the override keyword is not working!)
         */
        double g2g(const p_vec &p, Group &g1, Group &g2) {
          if (g1.id==Group::MOLECULAR) {
            if (g2.id==Group::MOLECULAR) {
              if ( geo->sqdist(g1.cm_trial, g2.cm_trial) > cut*cut ) {
                const GroupMolecular& m1 = static_cast<const GroupMolecular&>(g1);
                const GroupMolecular& m2 = static_cast<const GroupMolecular&>(g2);
                return 0; // v2v(m1.cg, m2.cg)
              }
            }
          }
          return Nonbonded<Tpairpot, Tgeometry>::g2g(p,g1,g2);
        }
      };

    /*!
     * \brief Energy class for hard-sphere overlap.
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

    //or simply add pointer to nonbonded<T>
    template<class Tpairpot, class Tgeometry>
      class Exclusions : public Nonbonded<Tpairpot,Tgeometry> {
        public:
          vector< std::pair<int,int> > pairs;
          Exclusions(InputMap &in) : Nonbonded<Tpairpot,Tgeometry>(in) {}
          //virtual double i2i(const p_vec &p, int i, int j) { return 0; }
      };

    /*!
     * \brief Class for handling bond pairs
     * \date Lund, 2011-2012
     * \author Mikael Lund
     *
     * Takes care of bonded interactions and can handle mixed bond types.
     * Example:
     * \code
     *    vector<particle> p(...);            // particle vector
     *    int i=10, j=11;                     // particle index
     *    Energy::Bonded b;
     *    b.add(i, j, Potential::Harmonic(0.1,5.0) );
     *    std::cout << b.info();
     *    double rij2 = ... ;                 // squared distance between i and j
     *    double u = b(i,j)( p[i], p[j], rij2 ); // i j bond energy in kT
     * \endcode
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
        double total(const p_vec&);
    };

    /*!
     * \brief Energy from external pressure for use in the NPT-ensemble.
     * \author Mikael Lund
     * \date Lund, 2011
     *
     * The system energy is:
     *
     * \f$\beta u = \beta pV - \ln V - N\ln V\f$.
     *
     * The two first terms are returned by external() while the last
     * term is obtained by summing g_external() over molecular
     * and or atomic groups.
     * If applied on an atomic group, \e N will be set to the number of
     * atoms in the group, while for a molecular group \e N =1.
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

    /*!
     * \brief External energy that will keep specific groups in a sub-volume of the system
     * \author Mikael Lund
     * \date Lund, 2012
     *
     * This energy class will check if particles in specific groups are located within a
     * rectangular box, spanned by two vector points, \c upper and \c lower. If outside
     * an infinite energy is returned. This is useful for constraining molecules in specific
     * parts of the simulation container.
     * Derived classes can re-implement the virtual
     * outside() function which should return \c true if a given point falls outside the
     * allowed region. Note that the current implementation can be problematic with
     * containers with periodic boundaries as the outside() function uses absolute positions.
     *
     * Example:
     * \code
     *   Energy::Hamiltonian pot;
     *   ...
     *   auto restricted = pot.create( Energy::RestrictedVolume(imap) );
     *   restricted->groups.push_back( &mygroup );
     * \endcode
     */
    class RestrictedVolume : public Energy::Energybase {
      private:
        Point upper, lower;
        string _info();
      protected:
        virtual bool outside(const Point&); //!< Determines if particle is outside allowed region
      public:
        std::vector<Group*> groups;              //!< List of groups to restrict
        RestrictedVolume(InputMap&, string="vconstrain");
        double g_external(const p_vec&, Group&) FOVERRIDE; //!< External energy working on group
    };

    /*!
     * \brief As Energy::RestrictedVolume but restrictions are applied only on the mass center
     * instead of all particles in group.
     */
    class RestrictedVolumeCM : public Energy::RestrictedVolume {
      public:
        RestrictedVolumeCM(InputMap&, string="vconstrain");
        double g_external(const p_vec&, Group&) FOVERRIDE; //!< External energy working on group
    };

    /*!
     * \brief Collection of Energybases that when summed give the Hamiltonian
     * \author Mikael Lund
     * \date Lund, 2011
     *
     * This class is used to collect several Energybase derivatives into a full hamiltonian.
     * The following example demonstrated how one can generate a Hamiltonian for bonded as
     * well as non-bonded interactions:
     * \code
     *   Energy::Hamiltonian pot;
     *   pot.create( Energy::Nonbonded<Tpairpot>(in) );
     *   pot.create( Energy::Bonded() );
     *   cout << pot.info();
     * \endcode
     *
     * Notice that we do not need to specify a Geometry for the Bonded energy class as this information
     * is simply passed on from the first added potential.
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

      /*!
       * \brief Create and add an energy class to energy list
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
    };

    /*!
     * \brief Constrain two group mass centra within a certain distance interval [mindist:maxdist]
     * \author Mikael Lund
     * \date Lund, 2012
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
        std::map< pair_permutable<Faunus::Group*>, data> gmap;
      public:
        void addPair(Group&, Group&, double, double);
        MassCenterConstrain(Geometry::Geometrybase&);
        double g_external(const p_vec&, Group&) FOVERRIDE;
    };

    /*!
     * \brief Dummy energy class that sums missed energy changes to avoid energy drifts
     * \author Mikael Lund
     *
     * This energy function is designed to be used with Move::Movebase classes that returns energy changes
     * not detectable in the energy drift checkup routines. The idea is simply to sum the energy change
     * discrepancy and treat this as an external potential. Use together with Energy::Hamiltonian.
     */
    class EnergyRest : public Energy::Energybase {
      private:
        double usum;
        string _info();
      public:
        EnergyRest();
        void add(double du); //!< Add energy change disrepancy, dU = U(metropolis) - U(as in drift calculation)
        double external() FOVERRIDE;
    };

    /*!
     * \brief Calculates the total system energy
     *
     * For a given particle vector, space, and energy class we try to calculate the
     * total energy taking into account inter- and intra-molecular interactions as well
     * as external potentials. While this may not work for all systems it may be a useful
     * first guess. This is the default energy routine for Move::ParallelTempering and may
     * also be used for checking energy drifts.
     */
    double systemEnergy(Space&, Energy::Energybase&, const p_vec&);

  }//Energy namespace
}//Faunus namespace
#endif
