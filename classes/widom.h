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
    unsigned long long int cnt;
    float u;
    particle a,b;
    container *con;
    average<float> expsum; 
    interaction<T_pairpot> *pot;
  public:
    widom(container &c, interaction<T_pairpot> &i,
        particle::type t1, particle::type t2) {
      cnt=0;
      con=&c;
      pot=&i;
      a=con->get(t1);
      b=con->get(t2);
    }
    string info();                              //!< Print results of analysis
    float muex() { return -log(expsum.avg()); } //!< Excess chemical potential
    float gamma() { return exp(muex()); }       //!< Activity coefficient
    
    //! Insert a salt pair and evaluate the excess chemical potential
    //! \param n Number of insertions
    void insert(unsigned char n=100) {
      for (unsigned char i=0; i<n; i++) {
        cnt++;
        con->randompos(a);
        con->randompos(b);
        u=pot->energy(con->p, a) +
          pot->energy(con->p, b) +
          pot->pair.pairpot(a,b)*pot->pair.f;
        expsum+=exp(-u);
      }
    }
};

template<class T>
string widom<T>::info() {
  ostringstream o;
  o << "# Widom Analysis:" << endl
    << "#   Number of insertions = " << cnt << endl
    << "#   Ion pair charges     = " << a.charge << ", " << b.charge << endl
    << "#   Excess chemical pot. = " << muex()  << endl
    << "#   Mean activity coeff. = " << gamma() << endl;
  return o.str();
}

#endif
