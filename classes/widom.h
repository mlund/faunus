#ifndef _WIDOM_H
#define _WIDOM_H

#include "average.h"
#include "container.h"
#include "potentials.h"

/*! \brief Widom method for excess chemical potentials
 *  \author Mikael Lund
 *  \todo Expand to one-particle insertion (Woodward+Svensson)
 *
 *  This class will insert "ghost" particles so as to
 *  calculate the mean excess chemical potential.
 */
template<class T_pairpot>
class widom {
  private:
    float u;
    particle a,b;
    container *con;
    average<float> expsum; 
    interaction<T_pairpot> *pot;
  public:
    widom(container &c, interaction<T_pairpot> &i,
        particle::type t1, particle::type t2) {
      con=&c;
      pot=&i;
      a=con->get(t1);
      b=con->get(t2);
    }
    float muex() { return -log(expsum.avg()); }
    float gamma() { return exp(muex()); }
    //! Insert ghost particles and evaluate excess chemical potential
    void insert(unsigned char n=100) {
      for (unsigned char i=0; i<n; i++) {
        con->randompos(a);
        con->randompos(b);
        u=pot->energy(con->p, a) +
          pot->energy(con->p, b) +
          pot->pair.pairpot(a,b)*pot->pair.f;
        expsum+=exp(-u);
      }
    }
};

#endif
