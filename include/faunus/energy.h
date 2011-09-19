#ifndef faunus_energy_base_h
#define faunus_energy_base_h

#include <iostream>
#include <vector>
#include <string>
#include <tuple>
#include <array>
#include <faunus/common.h>
#include <faunus/point.h>
#include <faunus/group.h>

// http://publib.boulder.ibm.com/infocenter/iadthelp/v8r0/index.jsp?topic=/com.ibm.xlcpp111.linux.doc/language_ref/variadic_templates.html
//
//
using namespace std;

namespace Faunus {

  namespace Energy {

    class energybase {
      protected:
        string name;
      public:
        static Geometry::geometrybase* geo; //!< Pointer to geometry functions - possibly slow, so avoid in time critical steps

        // single particle interactions
        virtual double all2all(const p_vec &p);
        virtual double i2i(const p_vec &p, int i, int j);
        virtual double i2g(const p_vec &p, const group &g, int i);
        virtual double i2all(const p_vec &p, int i);
        virtual double i_external(const p_vec &p, int i);
        virtual double i_internal(const p_vec &p, int i);

        // group interactions
        virtual double g2g(const p_vec &p, const group &g1, const group &g2);
        virtual double g2all(const p_vec &p, const group &g);
        virtual double g_external(const p_vec &p, const group &g);
        virtual double g_internal(const p_vec &p, const group &g);
        string info();
    };

    template<class Tpotential>
      class nonbonded : public energybase {
        private:
          typedef p_vec::const_iterator pviter;
        public:
          Tpotential pair;
          nonbonded(inputfile &in) : pair(in) {}

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
      };

    //or simply add pointer to nonbonded<T>
    template<class Tpotential>
      class exclusions : public nonbonded<Tpotential> {
        exclusions(inputfile &in) : nonbonded<Tpotential>(in) {}
        virtual double i2i(const p_vec &p, int i, int j) { return 0; }
      };

    template<class Tbondpot>
      class bonded : public energybase {
        //implement pointer to bond list
        // should we have a class that keeps track of particles,
        // groups, bond lists etc.?? YES!
        public:
          //void set_bondlist(int i, int j);
          double i_internal(const p_vec &p, int i) { return 0; }
          double g_internal(const p_vec &p, const group &g) { return 0; }
      };

    class hamiltonian : public energybase {
      private:
        vector<energybase*> base;
      public:
        hamiltonian(Geometry::geometrybase* geoPtr) { geo=geoPtr; }

        void operator+=(energybase &b) { base.push_back(&b); }

        // single particle interactions
        double i2i(const p_vec &p, int i, int j) {
          double u=0;
          for (int k=0; k<base.size(); ++k)
            u+=base[k]->i2i(p,i,j);
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
        string info() { return "hej";  };
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
  }//Energy namespace
}//Faunus namespace
#endif
