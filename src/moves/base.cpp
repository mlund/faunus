#include "faunus/moves/base.h"

#include "faunus/titrate.h"
#include "faunus/slump.h"
#include "faunus/io.h"

namespace Faunus {
  
  // interaction *pot;
  // pot=&inter; (pass "inter" by reference)
  markovmove::markovmove(ensemble &e, container &c, energybase &i) {
    du=utot=dp=deltadp=0;
    cnt=naccept=0;
    ens=&e;
    con=&c;
    pot=&i;
    runfraction=1;
  }

  string markovmove::info() {
    std::ostringstream o;
    o << endl
      << "# " << name << ":" << endl;
    if (cnt>0) {
      o << "#   Acceptance                = " << accepted()*100 << endl
        << "#   Number of trials          = " << cnt << endl
        << "#   Pct. of Markov steps      = " << runfraction*100 << endl
        << "#   Energy change (kT)        = " << utot << " " << utot/cnt << " "
                                        << utot/(accepted()*cnt) << endl;
      if (dp!=0) {
        o << "#   Displacement param.       = " << dp << endl;
        if (dpsqr.sum>0)
          o << "#   Average displacement      = " << sqrt(dpsqr.avg()) << endl
            << "#   Mean square displacement  = " << sqrt(dpsqr.sum) << endl;
      }
    }
    if (cite.empty()==false)
      o << "#   More information:           " << cite << endl;
    return o.str();
  }

  bool markovmove::run(float p) {
    return (p>slp.random_one()) ? true : false;
  }

  float markovmove::accepted() {
    return naccept/float(cnt);
  }

  /*!
   * This function will adjust the displacement parameter in a way
   * that the acceptance ration lies within a certain tange. Useful
   * for equilibration runs -- do not use it in production runs!
   * \param max Maximum percentage of accepted moves
   * \warning This violates the detailed balance criteria!
   * \param min Minimum percentage of accepted moves
   * \author Mikael Lund
   * \todo Specify a maxmimum dp
   */
  void markovmove::adjust_dp(float min, float max) {
    float a=accepted()*100.;
    if (a>max) dp+=deltadp;
    if (a<min) dp-=deltadp;
    if (dp<=0) dp=deltadp;
  }

}
