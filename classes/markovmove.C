#include "markovmove.h"

//-------------- SALT MOVE -------------
saltmove::saltmove(
    ensemble &e, container &c, interaction<T_pairpot> &i ) :
  markovmove(e,c,i) {}

/*! \param group Group containing mobile ions
 *  \param dp Displacement parameter
 *  \param n Particle number to move (-1 for random)
 */
void saltmove::move(group &g, float dp, int n) {
  cnt++;
  du=0;
  if (n==-1)
    n=g.random(); 
  con->displace(n, dp); 
  if (con->collision( con->trial[n] )==true)
    rc=HC;
  else {
    #pragma omp parallel
    {
      #pragma omp sections
      {
        #pragma omp section
        { uold = pot->energy(con->p, n);   }
        #pragma omp section
        { unew = pot->energy(con->trial,n);   }
      }
    }
    du = unew - uold;
    if (ens->metropolis(du)==true) {
      rc=OK;
      naccept++;
      con->p[n] = con->trial[n];
      return;
    }
    else
      rc=ENERGY;
    con->trial[n] = con->p[n];
  }
}
