#ifndef FAUNUS_BONDED_H
#define FAUNUS_BONDED_H

#include "faunus/energy/base.h"


// D R A F T ! ! !

namespace Faunus {

  class harmonicbond {
    public:
      unsigned int i;                //!< Particle 1
      unsigned int j;                //!< Particle 2
      double k;                      //!< Force constant
      double req;                    //!< Equilibrium distance
      double energy(container &);    //!< Perhaps a template function - needs container
      bool isMemberOf(unsigned int); //!< Test if particle is part of pair
  };
  bool harmonicbond::isMemberOf(unsigned int n) {
    return ( n==i || n==j ) ? true : false;
  }

  //--------------------------------------

  template <class Tbond> class bondlist {
    private:
      vector<Tbond> v;               //!< List of all bonds
      vector< vector<*Tbond> > bl;   //!< Bondlist for each particle in the system;

    public:
      void addBond(Tbond &b) {
        v.push_back(b);
      }

      void updateBondList(container &con) {
        bl.clear();
        bl.resize( con.p.size() );
        for (int b=0; b<v.size(); b++)
          for (int i=0; i<con.p.size(); i++)
            if (v[b].isMemberOf(i)==true)
              bl[i].push_back( &v[b] );
      }

      double energy(container &con, unsigned int i) {  //!< Energy of i'th particle
        double u=0;
        for (unsigned int b=0; b<bl[i].size(); b++)
          u+=bl[i][b]->energy(con);
        return u;
      }
  };

  //------------------------

  class bondEnergy {
    public:
      bondlist<harmonicbond> harmonic;
      bondlist<angularbond> angle;
  };

};//namespace
#endif
