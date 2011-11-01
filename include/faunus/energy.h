#ifndef faunus_energy_h
#define faunus_energy_h

#include <faunus/common.h>
#include <faunus/point.h>
#include <faunus/group.h>
#include <faunus/textio.h>
#include <faunus/potentials.h>

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
     * \brief General class for handling pairs of particles
     * \date Lund, October 2011
     *
     * This is a general class for handling properties for pairs. One example is bonds between
     * particles, identified through the two particle index. Properties are added through to
     * add() template function which can also handle derived classes of Tpairprop. Upon adding
     * new properties, space is dynamically allocated inside the class for each property object.
     */
    template<class Tpairprop, typename Tij=int>
      class ParticlePairs {
        typedef shared_ptr<Tpairprop> PropPtr;
        private:
        vector<PropPtr> created; //!< list of allocated pair properties
        protected:
        std::map<Tij, std::map<Tij, PropPtr> > list;
        string name;
        public:
        template<typename Tderived>
          void add(Tij i, Tij j, Tderived p) {
            assert(i!=j); //debug
            if (i!=j) {
              created.push_back( shared_ptr<Tderived>(new Tderived(p)) ); 
              list[i][j]=created.back();
              list[j][i]=created.back();
            }
          }

        //!< \brief Retrieve reference to bond via (i,j) operator
        Tpairprop& operator() (Tij i, Tij j) {
          assert( list[i][j]!=nullptr ); //debug
          return *list[i][j];
        }

        string info() {
          using namespace Faunus::textio;
          std::ostringstream o;
          o << indent(SUBSUB) << std::left
            << setw(7) << "i" << setw(7) << "j" << endl;
          for (auto i : list)
            for (auto j : list[i.first])
              if (i.first < j.first)
                o << indent(SUBSUB) << std::left << setw(7) << i.first
                  << setw(7) << j.first << j.second->brief() << endl;
          return o.str();
        }
      };

    /*!
     * \brief Class for handling bond pairs
     * \date Lund, October 2011
     * \author Mikael Lund
     *
     * Example:
     * \code
     *    vector<particle> p(...);            // particle vector
     *    int i=10, j=11;                     // particle index
     *    Energy::ParticleBonds bonds;
     *    bonds.add(i, j, Potential::Harmonic(0.1,5.0) );
     *    std::cout << bonds.info();
     *    double rij2 = ... ;                 // squared distance between i and j
     *    double u = bonds(i,j)( p[i], p[j], rij2 ); // i j bond energy
     * \endcode
     */
    class ParticleBonds : public ParticlePairs<Potential::PairPotentialBase,int> {
      typedef ParticlePairs<Potential::PairPotentialBase> pairs;
      public:
      ParticleBonds();
      double i2i(Geometry::Geometrybase&, const p_vec&, int, int);    //!< Bond energy of i'th particle with j'th
      double totalEnergy(Geometry::Geometrybase&, const p_vec&, int); //!< Bond energy of i'th particle (kT)
      double totalEnergy(Geometry::Geometrybase&, const p_vec&, const Group&); //!< Bond energy of group (kT)
      double totalEnergy(Geometry::Geometrybase&, const p_vec&); //!< Total bond energy of all bonds (kT)
    };

    /*!
     *  \brief Base class for energy evaluation
     *  \note All energy functions are expected to return energies in units of kT.
     *
     *  This base class defines functions for evaluating interactions between particles,
     *  groups, external potentials etc. By default all energy functions returns ZERO
     *  and derived classes are expected only to implement functions relevant for certain
     *  properties. I.e. a derived class for non-bonded interactions are not expected to
     *  implement i_internal(), for example.
     */
    class Energybase {
      private:
        virtual string _info()=0;
        char w; //!< Width of info output
      protected:
        Geometry::Geometrybase* geo; //!< Pointer to geometry used to calculate interactions
      public:
        string name;                                          //!< Short informative name
        Energybase();
        virtual ~Energybase();
        virtual Geometry::Geometrybase& getGeometry();        //!< Reference to geometry used for interactions
        bool setGeometry( Geometry::Geometrybase& );          //!< Set Geometrybase
        virtual double p2p(const particle&, const particle&); //!< Particle-particle energy
        virtual double all2p(const p_vec&, const particle&);  //!< Particle-Particle vector energy
        virtual double all2all(const p_vec&);                 //!< All inter-particle energies (N^2)
        virtual double i2i(const p_vec&, int, int);           //!< i'th particle with j'th particle
        virtual double i2g(const p_vec&, Group &, int);       //!< i'th particle with group
        virtual double i2all(const p_vec&, int);              //!< i'th particle with all other particles
        virtual double i_external(const p_vec&, int);         //!< internal energy of i'th particle
        virtual double i_internal(const p_vec&, int);         //!< External energy of i'th particle
        virtual double p_external(const particle&);           //!< External energy of particle
        double i_total(const p_vec&, int);                    //!< Total energy of i'th particle = i2all + i_external + i_internal
        virtual double g2g(const p_vec&, Group&, Group&);     //!< Group-Group energy
        virtual double g2all(const p_vec&, Group&);           //!< Energy of Group with all other particles
        virtual double g_external(const p_vec&, Group&);      //!< External energy of group
        virtual double g_internal(const p_vec&, Group&);      //!< Internal energy of group
        virtual double v2v(const p_vec&, const p_vec&);       //!< Particle vector-Particle vector energy
        virtual double external();                            //!< External energy - pressure, for example.
        string info();                                        //!< Information
    };

    /*!
     * \brief Energy class for non-bonded interactions.
     *
     * Tpotential is expected to be a pair potential with the following
     * properties:
     * \li pair(const InputMap&);
     * \li double pair.energy(const particle&, const particle&);
     * \li double pait.tokT();
     */
    template<class Tpotential>
      class Nonbonded : public Energybase {
        private:
          string _info() {
            return pair.info(25);
          }
        public:
          Tpotential pair;
          Nonbonded(InputMap &in) : pair(in) {
            name="Nonbonded N" + textio::squared + " - " + pair.name;
            geo=&pair.geo;
          }

          Geometry::Geometrybase& getGeometry() {
            geo=&pair.geo;
            return Energybase::getGeometry();
          }

          inline double p2p(const particle &a, const particle &b) {
            return pair.energy(a,b)*pair.tokT();
          }

          double all2p(const p_vec &p, const particle &a) {
            double u=0;
            for (auto &b : p)
              u+=pair.energy(a,b);
            return u*pair.tokT();
          }

          double all2all(const p_vec &p) {
            int n=p.size();
            double u=0;
            for (int i=0; i<n-1; ++i)
              for (int j=i+1; j<n; ++j)
                u+=pair.energy( p[i],p[j] );
            return u*pair.tokT();
          }

          double i2i(const p_vec &p, int i, int j) {
            return pair.energy(p[i],p[j]);
          }

          double i2g(const p_vec &p, Group &g, int j) {
            double u=0;
            int len=g.end+1;
            if (j>=g.beg && j<=g.end) {   //j is inside g - avoid self interaction
              for (int i=g.beg; i<j; i++)
                u+=pair.energy(p[i],p[j]);
              for (int i=j+1; i<len; i++)
                u+=pair.energy(p[i],p[j]);
            } else                        //simple - j not in g
              for (int i=g.beg; i<len; i++)
                u+=pair.energy(p[i],p[j]);
            return pair.tokT()*u;  
          }

          double i2all(const p_vec &p, int i) {
            double u=0;
            int n=(int)p.size();
            for (int j=0; j<i; ++j)
              u+=pair.energy( p[i], p[j] );
            for (int j=i+1; j<n; ++j)
              u+=pair.energy( p[i], p[j] );
            return u*pair.tokT();
          }

          double g2g(const p_vec &p, Group &g1, Group &g2) {
            double u=0;
            int ilen=g1.end+1, jlen=g2.end+1;
#pragma omp parallel for reduction (+:u) schedule (dynamic)
            for (int i=g1.beg; i<ilen; ++i)
              for (int j=g2.beg; j<jlen; ++j)
                u+=pair.energy(p[i],p[j]);
            return pair.tokT()*u;
          }

          double g2all(const p_vec &p, Group &g) {
            int ng=g.end+1, np=p.size();
            double u=0;
#pragma omp parallel for reduction (+:u)
            for (int i=g.beg; i<ng; ++i) {
              for (int j=0; j<g.beg; j++)
                u += pair.energy(p[i],p[j]);
              for (int j=ng; j<np; j++)
                u += pair.energy(p[i],p[j]);
            }
            return u*pair.tokT();
          }

          double g_internal(const p_vec &p, Group &g) { 
            if (g.beg==-1) return 0;
            double u=0;
            int step=1,n=g.end+1;
            for (int i=g.beg; i<n-step; i++)
              for (int j=g.beg+step*((i-g.beg)/step+1); j<n; j++)
                u+=pair.energy(p[i],p[j]);
            return pair.tokT()*u;
          }

      };

    template<class Tpotential>
      class Nonbonded_CG : public Nonbonded<Tpotential> {
        using Nonbonded<Tpotential>::geo;
        using Nonbonded<Tpotential>::name;
        public:
        double cut;
        Nonbonded_CG(InputMap &in) : Nonbonded<Tpotential>(in) {
          name+="(Molecular Group CG)";
        }

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
          return Nonbonded<Tpotential>::g2g(p,g1,g2);
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
          HardSphereOverlap(InputMap &mcp) : geo(mcp) {};
          inline double i2i(const p_vec &p, int i, int j) {
            return hs(p[i], p[j], geo.sqdist(p[i], p[j]) );
          }
          double all2all(const p_vec &p) {
            for (auto i=p.begin(); i!=p.end()-1; ++i)
              for (auto j=i+1; j!=p.end(); ++j)
                if ( hs(*i,*j, geo.sqdist(*i,&j) )==pc::infty )
                  return pc::infty;
            return 0;
          }
          double g2all(const p_vec &p, const Group &g) {
            for (int i=g.beg; i<=g.end; i++) {
              for (int j=0; j<g.beg; j++)
                if ( i2i(p,i,j)==pc::infty )
                  return pc::infty;
              for (int j=g.end+1; j<(int)p.size(); j++)
                if ( i2i(p,i,j)==pc::infty )
                  return pc::infty;
              return 0;
            }
          }
          double g2g(const p_vec &p, const Group &g1, const Group &g2) {
            for (int i=g1.beg; i<=g1.end; ++i)
              for (int j=g2.beg; j<=g2.end; ++j)
                if ( i2i(p,i,j)==pc::infty )
                  return pc::infty;
            return 0;
          }
      };

    //or simply add pointer to nonbonded<T>
    template<class Tpotential>
      class Exclusions : public Nonbonded<Tpotential> {
        public:
          vector< std::pair<int,int> > pairs;
          Exclusions(InputMap &in) : Nonbonded<Tpotential>(in) {}
          //virtual double i2i(const p_vec &p, int i, int j) { return 0; }
      };

    /*!
     * \brief Energy class for bonded interactions
     *
     * Takes care of bonded interactions and can handle mixed bond types.
     * Example:
     * \code
     * Energy::Bonded b(myGeometry);
     * Potential::Harmonic h(k, req);
     * b.bonds.add(10,12,h); // bond particle 10 and 12 
     * \endcode
     */
    class Bonded : public Energy::Energybase {
      private:
        string _info();
      public:
        Bonded();
        Bonded(Geometry::Geometrybase&);
        ParticleBonds bonds;
        //double i2i(const p_vec&, int, int);
        //double i2g(const p_vec&, Group&, int);
        double i2all(const p_vec&, int);
        double g_internal(const p_vec&, Group &);
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
     * term is obtained by summing g_external() over molecular groups.
     * If applied on an atomic group, N will be set to the number of
     * particles in the group.
     */
    class ExternalPressure : public Energy::Energybase {
      private:
        double P; //!< Pressure, p/kT
        string _info();
      public:
        ExternalPressure(Geometry::Geometrybase&, double);
        double external();  //!< External energy working on system. pV/kT-lnV
        double g_external(const p_vec&, Group&); //!< External energy working on group
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
      void setVolume(double);       //!< Set volume of all contained energy classes

      //!< \brief Create and add an energy class to energy list
      template<typename Tenergychild> shared_ptr<Tenergychild> create(Tenergychild c) {
        shared_ptr<Tenergychild> childptr( new Tenergychild(c) );
        childptr->getGeometry(); // not pretty...need to update geo pointer for i.e. nonbonded class
        created.push_back(childptr);
        add(*childptr);
        return childptr;
      }

      void add(Energybase&); //!< Add existing energy class to list
      double p2p(const particle&, const particle&);
      double all2p(const p_vec&, const particle&);
      double all2all(const p_vec&);
      double i2i(const p_vec&, int, int);
      double i2g(const p_vec&, Group&, int);
      double i2all(const p_vec&, int);
      double i_external(const p_vec&, int);
      double i_internal(const p_vec&, int);
      double g2g(const p_vec&, Group&, Group&);
      double g2all(const p_vec&, Group&);
      double g_external(const p_vec&, Group&);
      double g_internal(const p_vec&, Group&);
      double external();
    };

  }//Energy namespace
}//Faunus namespace
#endif
