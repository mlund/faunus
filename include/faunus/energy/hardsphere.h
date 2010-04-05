#ifndef FAUNUS_INT_HARDSPHERE_H
#define FAUNUS_INT_HARDSPHERE_H

#include "faunus/hardsphere.h"

namespace Faunus {
  /*!
   * \brief Hardsphere check, then normal potential function
   * \author Mikael Lund
   * \date Prague 2008
   * \todo Add overlap check in the system energy function.
   * \warning Untested!
   *
   * This interaction class first check for hardsphere overlap
   * and - if none found - proceeds with normal energy summation
   * according the specified pair potential.
   */
  template<class T>
  class interaction_hs : public interaction<T>, private hardsphere {
  private:
    double infty;
  public:
    interaction_hs(inputfile &in) : interaction<T>(in) { infty=1000.; }
    double energy(const particle &a, const particle &b) {
      return (a.overlap(b)==true) ? infty  : interaction<T>::energy(a,b);
    }
    double energy(vector<particle> &p) {
      return interaction<T>::energy(p);
    }
    double energy(vector<particle> &p, int i) {
      return (overlap(p,i)==true) ? infty : interaction<T>::energy(p,i);
    }
    double energy(vector<particle> &p, const particle &a) {
      return (overlap(p,a)==true) ? infty : interaction<T>::energy(p,a);
    }
    double energy(vector<particle> &p, const group &g) {
      return (overlap(p,g)==true) ? infty : interaction<T>::energy(p,g);
    }
  };
  
  template<class T>
  class interaction_vector : public interaction<T> {
  private:
    int len;
    double infty;
    double r2[4000], qq[4000];
  public:
    interaction_vector(inputfile &in) : interaction<T>(in) {
      interaction<T>::name+="Vectorized, Full N^2";
    }
    double energy(vector<particle> &p) {
      len=1;
      int n=p.size();
      for (int i=0; i<n-1; i++)
        for (int j=i+1; j<n; j++) {
          r2[len]   = interaction<T>::pair.sqdist(p[i],p[j]);
          qq[len++] = p[i].charge*p[j].charge;
        }
      return interaction<T>::pair.VectorEnergy(r2,qq,&len);
    }
    double energy(vector<particle> &p, const group &g) {
      len=1;
      int n=g.end+1, psize=p.size();
      for (int i=g.beg; i<n; ++i) {
        for (int j=0; j<g.beg; j++) {
          r2[len]   = interaction<T>::pair.sqdist(p[i],p[j]);
          qq[len++] = p[i].charge*p[j].charge;
        }
        for (int j=n; j<psize; j++) {
          r2[len]   = interaction<T>::pair.sqdist(p[i],p[j]);
          qq[len++] = p[i].charge*p[j].charge;
        }
      }
      return interaction<T>::pair.VectorEnergy(r2,qq,&len);
    }
    double energy(vector<particle> &p, const group &g1, const group &g2) {
      len=1;
      int ilen=g1.end+1, jlen=g2.end+1;
      //#pragma omp parallel for reduction (+:u)
      for (int i=g1.beg; i<ilen; i++)
        for (int j=g2.beg; j<jlen; j++) {
          r2[len]   = interaction<T>::pair.sqdist(p[i],p[j]);
          qq[len++]= p[i].charge*p[j].charge;
        }
      return interaction<T>::pair.VectorEnergy(r2,qq,&len);
    }
  };
  
}//namespace
#endif