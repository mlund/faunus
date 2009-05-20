#ifndef FAU_SPRINGENERGY_H
#define FAU_SPRINGENERGY_H

//#include "faunus/energy/base.h"

namespace Faunus {
  /*!
   * \brief Implementation of all energy functions
   * \author Mikael Lund
   *
   * This is an expression template that constructs the
   * energy functions of a particular pair potential. All
   * returned energies are in units of kT.
   *
   * \code
   * inputfile in("input.conf");
   * interaction<pot_coulomb> pot(in);
   * pot.energy(...);
   * \endcode
   */
  template<class T> class springinteraction : public interaction<T> {
    private:
      vector<unsigned short> l;
    public:
      using interaction<T>::pair;
      double k; //!< Spring constant
      double req; //!< Equilibrium distance
      springinteraction(inputfile &in) : interaction<T>(in) {
        interaction<T>::name="Full N^2 w. spring potential";
        k   = in.getflt("springconstant", 30);
        req = in.getflt("springeqdist", 2);
      };

      double u_monomer(const vector<particle> &p, const polymer &g, unsigned int i ) {
        // normal energy function if no neighbors
        if ( g.nb[i].size()==0 )
          return interaction<T>::energy(p,g,i);

        // energy w. all particles, skip neighbors (uses stl::find - maybe slow).
        double u=0;
        for (unsigned int j=0; j<i; j++)
          if ( g.areneighbors(i,j)==false )
            u+=pair.pairpot( p[i], p[j] );
        for (unsigned int j=i+1; j<p.size(); j++)
          if ( g.areneighbors(i,j)==false )
            u+=pair.pairpot( p[i], p[j] );
        u=u*pair.f;

        // calc. spring and coulomb energy w. neighbors
        for (unsigned int j=0; j<g.nb[i].size(); j++) {
          double r=sqrt( pair.sqdist( p[i], p[ g.nb[i][j] ] ) ),
                 dr = r - req;
          u+=k*dr*dr + pair.f * p[i].charge * p[g.nb[i][j]].charge / r;
        }
        return u;
      }
      double uself_polymer(const vector<particle> &p, const polymer &g) {
        double dr,r,u=0;
        int i,j,n=g.end+1;
        for (i=g.beg; i<n-1; i++)
          for (j=i+1; j<n; j++)
            if ( g.areneighbors(i,j)==false )
              u+=pair.f*pair.pairpot(p[i], p[j]);
            else {
              r=sqrt( pair.sqdist( p[i], p[j] ) );
              dr = r - req;
              u+=k*dr*dr + pair.f * p[i].charge*p[j].charge/r;
            }
        return u;
      }
  };
}//namespace
#endif
