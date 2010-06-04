#include "faunus/moves/rotational.h"

namespace Faunus {

  macrorot::macrorot( ensemble &e,
      container &c, energybase &i ) : markovmove(e,c,i)
  {
    name.assign("MACROMOLECULAR ROTATION");
    runfraction=1.0;
    deltadp=0.1;
    dp=1.0;
  };

  double macrorot::move(macromolecule &g) {
    if (slp.runtest(runfraction)==false)
      return 0;
    markovmove::move();
    g.rotate(*con, dp); 
    for (int i=g.beg; i<=g.end; i++) {
      if (con->collision(con->trial[i])==true) {
        g.undo(*con);
        dpsqr+=0;
        return 0;
      }
    }
    //#pragma omp parallel
    {
      //#pragma omp sections
      {
        //#pragma omp section
        { uold = pot->energy(con->p, g);   }
        //#pragma omp section
        { unew = pot->energy(con->trial, g);   }
      }
    }
    du = unew-uold;
    if (ens->metropolis(du)==true) {
      rc=OK;
      utot+=du;
      naccept++;
      g.accept(*con);
      return du;
    } else rc=ENERGY;
    du=0;
    g.undo(*con);
    return du;
  }

  double macrorot::move(vector<macromolecule> &g, int n) {
    if (slp.runtest(runfraction)==false)
      return 0;
    markovmove::move();
    g.at(n).rotate(*con, dp); 
    for (int i=g[n].beg; i<=g[n].end; i++) {
      if (con->collision(con->trial[i])==true) {
        g[n].undo(*con);
        dpsqr+=0;
        return 0;
      }
    }
    double deltau=0;
#pragma omp parallel for reduction (+:deltau) 
    for (int i=0; i<g.size(); i++)
      if (i!=n)
        deltau += pot->energy(con->trial, g[i], g[n]) - pot->energy(con->p, g[i], g[n]);
    du=deltau;
    if (ens->metropolis(du)==true) {
      rc=OK;
      utot+=du;
      naccept++;
      g[n].accept(*con);
      return du;
    } else rc=ENERGY;
    du=0;
    g[n].undo(*con);
    return du;
  }
}
