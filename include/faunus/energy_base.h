#ifndef faunus_energy_base_h
#define faunus_energy_base_h

#include <iostream>
#include <vector>
#include <string>
#include <tr1/tuple>
#include <tr1/array>
#include <faunus/point.h>

// http://publib.boulder.ibm.com/infocenter/iadthelp/v8r0/index.jsp?topic=/com.ibm.xlcpp111.linux.doc/language_ref/variadic_templates.html
//
//
using namespace std;

namespace Faunus {
  class group;

  namespace Energy {

    class energybase {
      public:
        // single particle interactions
        virtual double i2i(const vector<particle> &p, int i, int j) { return 0; }
        virtual double i2g(const vector<particle> &p, const group &g, int i) { return 0; }
        virtual double i2all(const vector<particle> &p, int i) { return 0; }
        virtual double i_external(const vector<particle> &p, int i) { return 0; }
        virtual double i_internal(const vector<particle> &p, int i) { return 0; }

        // group interactions
        virtual double g2g(const vector<particle> &p, const group &g1, const group &g2) { return 0; }
        virtual double g2all(const vector<particle> &p, const group &g) { return 0; }
        virtual double g_external(const vector<particle> &p, const group &g) { return 0; }
        virtual double g_internal(const vector<particle> &p, const group &g) { return 0; }
        string info() { return "hej";  };
    };

    template<class T>
      class nonbonded : public energybase {
        public:
          T pair;
          nonbonded(inputfile &in) : pair(in) {
          }
          virtual double i2i(const vector<particle> &p, int i, int j) { return 0; }
          virtual double i2g(const vector<particle> &p, const group &g, int i) { return 0; }
          virtual double i2all(const vector<particle> &p, int i) { return 0; }
          virtual double i_internal(const vector<particle> &p, int i) { return 0; }

          virtual double g2g(const vector<particle> &p, const group &g1, const group &g2) { return 0; }
          virtual double g2all(const vector<particle> &p, const group &g) { return 0; }
          virtual double g_internal(const vector<particle> &p, const group &g) { return 0; }
      };

    class bonded : public energybase {
      //implement pointer to bond list
      // should we have a class that keeps track of particles,
      // groups, bond lists etc.?? YES!
      public:
        //void set_bondlist(int i, int j);
        double i_internal(const vector<particle> &p, int i) { return 0; }
        double g_internal(const vector<particle> &p, const group &g) { return 0; }
    };

    class hamiltonian : public energybase {
      private:
        vector<energybase*> base;
      public:
        void operator+=(energybase &b) {
          base.push_back(&b);
        }

        // single particle interactions
        double i2i(const vector<particle> &p, int i, int j) {
          double u=0;
          for (int k=0; k<base.size(); ++k)
            u+=base[k]->i2i(p,i,j);
          return u;
        }

        double i2g(const vector<particle> &p, const group &g, int i) { return 0; }
        double i2all(const vector<particle> &p, int i) { return 0; }
        double i_external(const vector<particle> &p, int i) { return 0; }
        double i_internal(const vector<particle> &p, int i) { return 0; }

        // group interactions
        double g2g(const vector<particle> &p, const group &g1, const group &g2) { return 0; }
        double g2all(const vector<particle> &p, const group &g) { return 0; }
        double g_external(const vector<particle> &p, const group &g) { return 0; }
        double g_internal(const vector<particle> &p, const group &g) { return 0; }
        string info() { return "hej";  };
    };
  }//namespace
}//namespace
#endif
