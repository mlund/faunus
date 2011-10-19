#ifndef faunus_energy_h
#define faunus_energy_h

#include <faunus/common.h>
#include <faunus/point.h>
#include <faunus/group.h>
#include <faunus/textio.h>
#include <faunus/bonded.h>

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
     *  groups, external potentials etc. By default all energy functions returns ZERO
     *  and derived classes are expected only to implement functions relevant for certain
     *  properties. I.e. a derived class for non-bonded interactions are not expected to
     *  implement i_internal(), for example.
     */
    class Energybase {
      protected:
        Geometry::Geometrybase* geo; //!< Pointer to geometry used to calculate interactions
      public:
        string name;

        Energybase();
        virtual ~Energybase();
        virtual Geometry::Geometrybase& getGeometry(); //!< Return reference to Geometrybase used for interactions
        bool setGeometry( Geometry::Geometrybase& ); //!< Set Geometrybase

        // interaction with external particle
        virtual double p2p(const particle&, const particle&);
        virtual double all2p(const p_vec&, const particle&);

        // single particle interactions
        virtual double all2all(const p_vec&);
        virtual double i2i(const p_vec&, int, int);
        virtual double i2g(const p_vec&, Group &, int);
        virtual double i2all(const p_vec&, int);
        virtual double i_external(const p_vec&, int); //!< Internal energy of i'th particle
        virtual double i_internal(const p_vec&, int);
        virtual double p_external(const particle&);
        double i_total(const p_vec&, int); //!< Total energy = i2all + i_external + i_internal

        // Group interactions
        virtual double g2g(const p_vec&, Group&, Group&);
        virtual double g2all(const p_vec&, Group&);
        virtual double g_external(const p_vec&, Group&);
        virtual double g_internal(const p_vec&, Group&);

        virtual double v2v(const p_vec&, const p_vec&);

        virtual string info();
     };

    template<class Tpotential>
      class Nonbonded : public Energybase {
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
            return pair.energy(a,b)*pair.tokT;
          }

          double all2p(const p_vec &p, const particle &a) {
            double u=0;
            for (auto &b : p)
              u+=pair.energy(a,b);
            return u*pair.tokT;
          }

          double all2all(const p_vec &p) {
            int n=p.size();
            double u=0;
            for (int i=0; i<n-1; ++i)
              for (int j=i+1; j<n; ++j)
                u+=pair.energy( p[i],p[j] );
            return u*pair.tokT;
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
            return pair.tokT*u;  
          }

          double i2all(const p_vec &p, int i) {
            double u=0;
            int n=(int)p.size();
            for (int j=0; j<i; ++j)
              u+=pair.energy( p[i], p[j] );
            for (int j=i+1; j<n; ++j)
              u+=pair.energy( p[i], p[j] );
            return u*pair.tokT;
          }

          double g2g(const p_vec &p, Group &g1, Group &g2) {
            double u=0;
            int ilen=g1.end+1, jlen=g2.end+1;
#pragma omp parallel for reduction (+:u) schedule (dynamic)
            for (int i=g1.beg; i<ilen; ++i)
              for (int j=g2.beg; j<jlen; ++j)
                u+=pair.energy(p[i],p[j]);
            return pair.tokT*u;
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
            return u*pair.tokT;
          }

          double g_internal(const p_vec &p, Group &g) { 
            if (g.beg==-1) return 0;
            double u=0;
            int step=1,n=g.end+1;
            for (int i=g.beg; i<n-step; i++)
              for (int j=g.beg+step*((i-g.beg)/step+1); j<n; j++)
                u+=pair.energy(p[i],p[j]);
            return pair.tokT*u;
          }

          string info() {
            std::ostringstream o;
            o << Energybase::info() << pair.info(25);
            return o.str();
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
                const Molecular& m1 = static_cast<const Molecular&>(g1);
                const Molecular& m2 = static_cast<const Molecular&>(g2);
                return 0; // v2v(m1.cg, m2.cg)
              }
            }
          }
          return Nonbonded<Tpotential>::g2g(p,g1,g2);
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

    class Bonded : public Energy::Energybase {
      public:
        Bonded();
        Bonded(Geometry::Geometrybase&);
        Faunus::Energy::ParticleBonds bonds;
        double i_internal(const p_vec &p, int i);
        double g_internal(const p_vec &p, Group &g);
        string info();
    };

    /*!
     * \brief Collection of Energybases that when summed gives the Hamiltonian
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
     * \endcode
     */
    class Hamiltonian : public Energybase {
      typedef shared_ptr<Energybase> baseptr;
      private:
      vector<baseptr> created;      //!< smart pointer list of created energy children
      public:
      vector<Energybase*> baselist; //!< pointer list to energy children to be summed

      template<typename Tenergychild> shared_ptr<Tenergychild> create(Tenergychild c) {
        shared_ptr<Tenergychild> childptr( new Tenergychild(c) );
        childptr->getGeometry(); // not pretty...need to update geo pointer for i.e. nonbonded class
        created.push_back(childptr);
        add(*childptr);
        return childptr;
      }

      void add(Energybase&); //!< Add existing energybase to list

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
      string info();
    };

  }//Energy namespace
}//Faunus namespace
#endif
