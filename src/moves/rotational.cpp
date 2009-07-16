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

  /*!
   * \todo Cell overlap test missing
   */
  double macrorot::move(macromolecule &g) {
    du=0;
    if (slp.runtest(runfraction)==false)
      return du;
    cnt++;
    g.rotate(*con, dp); 
    //insert cell overlap test
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

}
