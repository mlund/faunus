#include <faunus/move_base.h>
#include <faunus/slump.h>

namespace Faunus {
  
  bool mcmove::metropolis(double du) {
    if (du > 0)
      if ( slp.random_one()>exp(-du) )
        return false;
    return true;
  }
  
  void mcmove::trialmove() {
    cnt++;
  }

  void mcmove::acceptmove() {
    cnt_accepted++;
  }
  
  double mcmove::energychange() {
    return 0;
  }
  
  double mcmove::move() {
    if (runfraction>slp.random_one())
      return 0;
    trialmove();
    double du=energychange();
    if (metropolis(du)==true) {
      acceptmove();
      dusum+=du;
    } else {
      rejectmove();
      du=0;
    }
    return du;
  }
  
}//namespace