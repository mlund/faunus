#ifndef _ANALYSIS_H
#define _ANALYSIS_H

#include "average.h"
#include "container.h"
#include "potentials.h"

class analysis {
  public:
    virtual string info()=0;
};

class systemenergy : public analysis {
  private:
    double u0,sum,cur;
    average<double> uavg, u2avg;
  public:
    string info();
    systemenergy(double);
    void setcurrent(double);    //!< Tell about the current system energy
    void operator+=(double);    //!< Add system energy change
};

string systemenergy::info() {
  uavg+=sum;
  u2avg+=sum*sum;
  ostringstream o;
  o << "# System energy (kT):" << endl
    << "#   Averages <U> <U^2> = " << uavg.avg() << " " << u2avg.avg() << endl
    << "#   Initial energy     = " << u0 << endl
    << "#   Initial + changes  = " << sum << endl
    << "#   Current energy     = " << cur << endl
    << "#   Absolute drift     = " << abs(cur-sum) << endl;
  return o.str();
}

/*! \brief Widom method for excess chemical potentials
 *  \author Mikael Lund
 *  \todo Expand to one-particle insertion (Woodward+Svensson)
 *
 *  This class will insert "ghost" particles so as to
 *  calculate the mean excess chemical potential.
 */
template<class T_pairpot>
class widom : public analysis {
  private:
    unsigned long long int cnt;
    double u;
    particle a,b;
    container *con;
    average<double> expsum; 
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
    double muex() { return -log(expsum.avg()); }//!< Excess chemical potential
    double gamma() { return exp(muex()); }      //!< Activity coefficient
    void insert(unsigned char=100);             //!< Widom insertions
};

//! Insert a salt pair and evaluate the excess chemical potential
//! \param n Number of insertions
template<class T>
void widom<T>::insert(unsigned char n) {
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
