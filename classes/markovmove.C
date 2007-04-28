#include "markovmove.h"

//-------------- SALT MOVE -------------
saltmove::saltmove(
    ensemble &e, space &s,
    interaction<T_pairpot> &i, container &c ) : markovmove(e,s,i)
{
  cPtr=&c;
}

/*! \param group Group containing mobile ions
 *  \param dp Displacement parameter
 *  \param n Particle number to move (-1 for random)
 */
void saltmove::move(group &g, float dp, int n) {
  cnt++;
  du=0;
  if (n>-1)
    n=g.random(); 
  s->displace(n, dp); 
  if (cPtr->collision( s->trial[n] )==true)
    rc=HC;
  else {
    #pragma omp parallel
    {
      #pragma omp sections
      {
        #pragma omp section
        { uold = pot->energy(s->p, n);   }
        #pragma omp section
        { unew = pot->energy(s->trial,n);   }
      }
    }
    du = unew - uold;
    if (ensPtr->metropolis(du)==true) {
      rc=OK;
      naccept++;
      s->p[n] = s->trial[n];
    }
    else
      rc=ENERGY;
  }
}
