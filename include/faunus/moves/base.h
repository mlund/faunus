#ifndef FAU_markovmove_h
#define FAU_markovmove_h

#include "faunus/common.h"
#include "faunus/container.h"
#include "faunus/energy.h"
#include "faunus/ensemble.h"
#include "faunus/titrate.h"
#include "faunus/slump.h"
#include "faunus/io.h"
#include "faunus/histogram.h"

namespace Faunus {
  /*! \brief Base class for MC moves
   *  \author Mikael Lund
   *
   *  This class will keep track of the returned energy, if the
   *  move was successful of not and can provide statistics about the
   *  accepted and rejected moves.
   *
   *  Each derived class should provide a simple move() function that will
   *  perform a trial move, calculate the energy and either accept or reject
   *  the move. Unsuccessful moves are automatically undone.
   */
  class markovmove {
    protected:
      slump slp;
      float runfraction;
      double uold, unew, deltadp;
      unsigned long long int cnt, naccept;
      string name;                        //!< Arbitrary name for the move
      string cite;                        //!< Reference to additional info. (article, url etc.)
      container *con;
      ensemble *ens;
      average<float> dpsqr;               //!< Average displacement squared
    public:
      enum keys {OK, ENERGY, HC};
      keys rc;                            //!< Return code from move() functions
      double dp,                          //!< Displacement parameter
             du,                          //!< Energy change of last move
             utot;                        //!< Sum of energy changes for all moves
      float accepted();                   //!< Return fraction of accepted moves
      bool run(float);                    //!< Probability
      void adjust_dp(float=30, float=40); //!< Adjust displacement parameter
      virtual string info();              //!< Show info about group 
      energybase *pot;
      markovmove(ensemble &, container &, energybase &);
  };

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
      o << "#   Acceptance          = " << accepted()*100 << endl
        << "#   Number of trials    = " << cnt << endl
        << "#   Pct. of Markov steps= " << runfraction*100 << endl
        << "#   Total energy change = " << utot << endl;
      if (dp!=0) {
        o << "#   Displacement param. = " << dp << endl
          << "#   Average displacement= " << sqrt(dpsqr.avg()) << endl;
      }
    }
    if (cite.empty()==false)
      o << "#   More information:     " << cite << endl;
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
#endif
