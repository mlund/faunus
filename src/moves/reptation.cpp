#include "faunus/moves/reptation.h"
#include "faunus/energy/base.h"
#include "faunus/ensemble.h"

namespace Faunus {

  reptation::reptation(ensemble &e, container &c, energybase &i, inputfile &in) : markovmove(e,c,i) {
    name.assign("Reptation move");
    prefix.assign("reptation_");
    runfraction=1.0;
    deltadp=0.;
  }

  double reptation::move(polymer &g) {  //This move assumes all bonds are equidistant!
    du=0;
    cnt++;
    point uv;
    uv.ranunit(slp);
    double bd;
    bd=con->dist(con->trial[g.beg],con->trial[g.beg+1]);

    if (slp.random_one()<0.5) {//Move forward
      direction=FORWARD;
      for(int i=g.beg; i<g.end; i++) {
        con->trial[i].x=con->trial[i+1].x;
        con->trial[i].y=con->trial[i+1].y;
        con->trial[i].z=con->trial[i+1].z;
      }
      con->trial[g.end]+=uv*bd;
    } else {                   //Move backward
      direction=BACKWARD;
      for(int i=g.end; i>g.beg; i--) {
        con->trial[i].x=con->trial[i-1].x;
        con->trial[i].y=con->trial[i-1].y;
        con->trial[i].z=con->trial[i-1].z;
      }
      con->trial[g.beg]+=uv*bd;
    }

    bool hc=false;
    if (con->collision(con->trial[g.beg])==true || con->collision(con->trial[g.end])==true) {
      hc=true;
      } else {
         unew=pot->energy(con->trial,g) + pot->internal(con->trial,g);
         uold=pot->energy(con->p,g) + pot->internal(con->p,g);
         du =  unew-uold;
 
      }
    if (hc==true) {
      rc=HC;
      du=0;
      g.undo(*con);
      return du;
    }
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

  string reptation::info() {
    std::ostringstream o;
    if (runfraction>0) {
      o << markovmove::info();
    }
    return o.str();
  }
} // namespace
