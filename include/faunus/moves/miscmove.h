#ifndef FAU_ROTTRANSLATIONAL_H
#define FAU_ROTTRANSLATIONAL_H
#include "faunus/moves/base.h"
#include "faunus/moves/translational.h"
namespace Faunus {

  /*! \brief Rotate group around cm plus dr (effectively a combined rotation and translation).
   *  \author Bjorn Persson
   *  \date Lund 2008
   */
  class transrot : public markovmove { 
    public:
      transrot( ensemble&, container&, energybase&);
      double move(macromolecule &);
      double dpt, dpr;      // step parameters
      string info();
      
  };
  
  transrot::transrot( ensemble &e,
      container &c, energybase &i ) : markovmove(e,c,i)
  {
    name = "COMBINED ROTATION and TRANSLATION";
    runfraction=1.0;
    deltadp=0.1;
    dp=1.0;
    dpt=dpr=0.5;  
  };

  string transrot::info() {
    std::ostringstream o;
    o <<  markovmove::info()
      << "#   Dp (translation)      = " <<  dpt << endl
      << "#   Dp (rotation)         = " <<  dpr << endl;
    return o.str();
  } 
  /*!
   * \todo Cell overlap test missing, dp is used both for dr and the angle
   */
  double transrot::move(macromolecule &g) {
    du=0;
    cnt++;
    g.rotate(*con, dpr, dpt); //dpt in ang, dpr in rad./2
    //insert cell overlap test
    for (int i=g.beg; i<=g.end; i++) { 
      if (con->collision( con->trial[i] )==true) { 
        rc=HC;
        dpsqr+=0.;
        g.undo(*con);
        return du; 
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
      dpsqr+=con->sqdist(g.cm, g.cm_trial);
      g.accept(*con);
      return du;
    } else rc=ENERGY;
    du=0;
    dpsqr+=0.;
    g.undo(*con);
    return du;
  }
}
#endif
