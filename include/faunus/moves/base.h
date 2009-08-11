#ifndef FAU_markovmove_h
#define FAU_markovmove_h

#include "faunus/common.h"
#include "faunus/container.h"
#include "faunus/energy.h"
#include "faunus/ensemble.h"

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
      double uold, unew, deltadp;
      unsigned long long int cnt, naccept;
      string cite;                        //!< Reference to additional info. (article, url etc.)
      container *con;
      ensemble *ens;
      average<float> dpsqr;               //!< Average displacement squared
    public:
      string name;                        //!< Arbitrary name for the move
      enum keys {OK, ENERGY, HC};
      keys rc;                            //!< Return code from move() functions
      double dp,                          //!< Displacement parameter
             du,                          //!< Energy change of last move
             utot;                        //!< Sum of energy changes for all moves
      float runfraction;                  //!< Fractional chance that the move will be performed
      float accepted();                   //!< Return fraction of accepted moves
      bool run(float);                    //!< Probability
      void adjust_dp(float=30, float=40); //!< Adjust displacement parameter
      void check(checkValue &);           //!< Output test
      virtual string info();              //!< Show info about group 
      energybase *pot;
      markovmove(ensemble &, container &, energybase &);
  };
}//namespace
#endif
