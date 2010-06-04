#ifndef FAU_SPRINGENERGY_H
#define FAU_SPRINGENERGY_H

#include "faunus/energy/base.h"

namespace Faunus {
  /*!
   * \brief Implementation of energy functions bonded interactions between
   *        polymers
   * \author Mikael Lund
   *
   * This template included interactions for harmonic bonds in
   * bonded polymers. The potential class given must have a function
   * double bond(particle&, particle&) that calculates the energy
   * between two neighboing atoms.
   */
  template<class T> class springinteraction : public interaction<T> {
    private:
      vector<unsigned short> l;
    public:
      using interaction<T>::pair;
      springinteraction(inputfile &in) : interaction<T>(in) {
        interaction<T>::name="Full N^2 w. spring potential";
      }

      /*!
       * Bonded and non-bonded energy of monomer i with the rest of the system
       */
      double u_monomer(vector<particle> &p, const polymer &g, unsigned int i ) {
        // normal energy function if no neighbors
        if ( g.nb[i-g.beg].size()==0 )
          return interaction<T>::energy(p,i);

        // non-bonded energy - skip neighbors
        double u=0;
        for (int j=0; j<i; j++) 
          if ( g.areneighbors(i,j)==false ) 
            u+=pair.pairpot( p[i], p[j] );
        for (int j=i+1; j<p.size(); j++) 
          if ( g.areneighbors(i,j)==false )
            u+=pair.pairpot( p[i], p[j] );
        u=u*pair.f;

        // bonded energy w. neighbors
        // - decrease non-bonded repulsion
        for (int j=0; j<g.nb[i-g.beg].size(); j++)
          u+=pair.f*pair.bond(p[i],p[ g.nb[i-g.beg][j]]);
        return u;
      }

      double uself_polymer(vector<particle> &p, const polymer &g) {
        double u=0;
        int n=g.end+1;
        for (int i=g.beg; i<n-1; i++)
          for (int j=i+1; j<n; j++)
            if ( g.areneighbors(i,j)==false )
              u+=pair.f*pair.pairpot(p[i], p[j]);
            else
              u+=pair.f*pair.bond(p[i],p[j]);
        return u;
      }

      double uself_popscmem(vector<particle> &p, const popscmembrane &g) {
        double u=0;
        int i,j;
        for (i=0; i<g.pops.size(); i++) {
          u+=uself_polymer(p, g.pops[i]);
          for (j=i+1; j<g.pops.size(); j++) 
            u+=interaction<T>::energy(p,g.pops[i],g.pops[j]);
          for (j=0; j<g.popc.size(); j++)
            u+=interaction<T>::energy(p,g.pops[i],g.popc[j]);
        }
        for (i=0; i<g.popc.size(); i++) {
          u+=uself_polymer(p, g.popc[i]); 
          for (j=i+1; j<g.popc.size();j++)
            u+=interaction<T>::energy(p, g.popc[i], g.popc[j]);
        }
        return u;
      }
  };
}//namespace
#endif
