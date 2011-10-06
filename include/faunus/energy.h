#ifndef faunus_energy_base_h
#define faunus_energy_base_h

#include <faunus/common.h>
#include <faunus/point.h>
#include <faunus/group.h>
#include <faunus/species.h>
#include <faunus/titrate.h>
#include <faunus/textio.h>

// http://publib.boulder.ibm.com/infocenter/iadthelp/v8r0/index.jsp?topic=/com.ibm.xlcpp111.linux.doc/language_ref/variadic_templates.html
//
//
using namespace std;

namespace Faunus {

  namespace Energy {

    class Energybase {
      protected:
        Geometry::Geometrybase* geo; //!< Pointer to geometry used to calculate interactions
      public:
        string name;

        Energybase();
        virtual ~Energybase();

        Geometry::Geometrybase& getGeometry(); //!< Return reference to Geometrybase used for interactions

        // interaction with external particle
        virtual double p2p(const particle&, const particle&);
        virtual double all2p(const p_vec&, const particle&);

        // single particle interactions
        virtual double all2all(const p_vec&);
        virtual double i2i(const p_vec&, int, int);
        virtual double i2g(const p_vec&, const group &, int);
        virtual double i2all(const p_vec&, int);
        virtual double i_external(const p_vec&, int);
        virtual double i_internal(const p_vec&, int);

        // group interactions
        virtual double g2g(const p_vec&, const group&, const group&);
        virtual double g2all(const p_vec&, const group&);
        virtual double g_external(const p_vec&, const group&);
        virtual double g_internal(const p_vec&, const group&);

        string info();
    };

    template<class Tpotential>
      class Nonbonded : public Energybase {
        private:
          typedef p_vec::const_iterator pviter;
        public:
          Tpotential pair;
          Nonbonded(InputMap &in) : pair(in) {
            name="Nonbonded N" + textio::squared + " - " + pair.name;
            geo=&pair.geo;
          }

          inline double p2p(const particle &a, const particle &b) {
            return pair.pairpot(a,b)*pair.tokT;
          }

          double all2p(const p_vec &p, const particle &a) {
            double u=0;
            for (auto &b : p)
              u+=pair.pairpot(a,b);
            return u*pair.tokT;
          }

          double all2all(const p_vec &p) {
            int n=p.size();
            double u=0;
            for (int i=0; i<n-1; ++i)
              for (int j=i+1; j<n; ++j)
                u+=pair.pairpot( p[i],p[j] );
            return u*pair.tokT;
          }

          double i2i(const p_vec &p, int i, int j) {
            return pair.pairpot(p[i],p[j]);
          }

          double i2g(const p_vec &p, const group &g, int j) {
            double u=0;
            int len=g.end+1;
            if (j>=g.beg && j<=g.end) {   //j is inside g - avoid self interaction
              for (int i=g.beg; i<j; i++)
                u+=pair.pairpot(p[i],p[j]);
              for (int i=j+1; i<len; i++)
                u+=pair.pairpot(p[i],p[j]);
            } else                        //simple - j not in g
              for (int i=g.beg; i<len; i++)
                u+=pair.pairpot(p[i],p[j]);
            return pair.tokT*u;  
          }

          double i2all(const p_vec &p, int i) {
            double u=0;
            for (int j=0; j<i; ++j)
              u+=pair.pairpot( p[i], p[j] );
            for (int j=i+1; j<p.size(); ++j)
              u+=pair.pairpot( p[i], p[j] );
            return u*pair.tokT;
          }
          double i_internal(const p_vec &p, int i) { return 0; }

          double g2g(const p_vec &p, const group &g1, const group &g2) {
            double u=0;
            int ilen=g1.end+1, jlen=g2.end+1;
#pragma omp parallel for reduction (+:u) schedule (dynamic)
            for (int i=g1.beg; i<ilen; ++i)
              for (int j=g2.beg; j<jlen; ++j)
                u+=pair.pairpot(p[i],p[j]);
            return pair.tokT*u;
          }

          double g2all(const p_vec &p, const group &g) {
            int ng=g.end+1, np=p.size();
            double u=0;
#pragma omp parallel for reduction (+:u)
            for (int i=g.beg; i<ng; ++i) {
              for (int j=0; j<g.beg; j++)
                u += pair.pairpot(p[i],p[j]);
              for (int j=ng; j<np; j++)
                u += pair.pairpot(p[i],p[j]);
            }
            return u*pair.tokT;
          }

          double g_internal(const p_vec &p, const group &g) { 
            if (g.beg==-1) return 0;
            double u=0;
            int step=1,n=g.end+1;
            for (int i=g.beg; i<n-step; i++)
              for (int j=g.beg+step*((i-g.beg)/step+1); j<n; j++)
                u+=pair.pairpot(p[i],p[j]);
            return pair.tokT*u;
          }

          string info() {
            using namespace textio;
            std::ostringstream o;
            o << header(name) << pair.info(25);
            return o.str();
          }
      };

    //or simply add pointer to nonbonded<T>
    template<class Tpotential>
      class Exclusions : public Nonbonded<Tpotential> {
        Exclusions(InputMap &in) : Nonbonded<Tpotential>(in) {}
        virtual double i2i(const p_vec &p, int i, int j) { return 0; }
      };

    template<class Tbondpot>
      class Bonded : public Energybase {
        //implement pointer to bond list
        // should we have a class that keeps track of particles,
        // groups, bond lists etc.?? YES!
        public:
          //void set_bondlist(int i, int j);
          double i_internal(const p_vec &p, int i) { return 0; }
          double g_internal(const p_vec &p, const group &g) { return 0; }
      };

    /*!
     * \brief Collection of Energybases that when summed gives the Hamiltonian
     * \author Mikael Lund
     * \date Lund, 2011
     */
    class Hamiltonian : public Energybase {
      protected:
        vector<Energybase*> baselist;

      public:

        void add(Energybase &e) {
          assert(&e!=NULL);
          geo=&e.getGeometry();
          baselist.push_back(&e);
        }

        // single particle interactions
        double i2i(const p_vec &p, int i, int j) {
          double u=0;
          for (auto b : baselist)
            u += b->i2i( p,i,j );
          return u;
        }

        double i2g(const p_vec &p, const group &g, int i) { return 0; }
        double i2all(const p_vec &p, int i) { return 0; }
        double i_external(const p_vec &p, int i) { return 0; }
        double i_internal(const p_vec &p, int i) { return 0; }

        // group interactions
        double g2g(const p_vec &p, const group &g1, const group &g2) { return 0; }
        double g2all(const p_vec &p, const group &g) { return 0; }
        double g_external(const p_vec &p, const group &g) { return 0; }
        double g_internal(const p_vec &p, const group &g) { return 0; }

        string info() {
          using namespace textio;
          char w=25;
          std::ostringstream o;
          o << header("Hamiltonian");
          o << pad(SUB,w,"Number of energybases") << baselist.size() << endl;
          for (auto e : baselist)
            o << pad(SUB,w,e->name) << endl;
          return o.str();
        }
    };

    class bondlist {
      public:
        enum bondtype {HARMONIC};
        class bond {
          public:
            int first,second;
            double req;
            double k;
            bool contain(const int &i) const { return (i==first||i==second) ? true:false; } 
        };
        bond& operator[](const int &i) {return v.at(i);}
        int size() const { return v.size(); }
        void add(int first, int second, double req, double k) {
          if (first>second)
            std::swap(second,first);
          bond b={first,second,req,k};
          v.push_back(b);
        }
      private:
        vector<bond> v;
    };

    class eqenergy : public Energybase {
      protected:
        typedef std::map<short, double> Tmap;
        typedef Tmap::iterator mapiter;
        Tmap energymap;                //!< Map of intrinsic energy for titratable sites
      public:
        eqenergy(equilibria &eq) {
          for (auto &s : atom.list )
            for (auto &process : eq.process )
              if ( process.one_of_us( s.id ) )
                energymap[s.id] = process.energy( s.id );
        }
        double i_internal(const p_vec &p, int i) {
          if (energymap.find( p[i].id )!=energymap.end())
            return energymap[ p[i].id ];
          return 0;
        }
    };

  }//Energy namespace
}//Faunus namespace
#endif
