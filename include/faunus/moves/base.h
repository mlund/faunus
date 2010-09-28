#ifndef FAU_markovmove_h
#define FAU_markovmove_h

#include "faunus/common.h"
#include "faunus/xytable.h"
#include "faunus/titrate.h"
#include "faunus/slump.h"

namespace Faunus {
  class container;
  class ensemble;
  class energybase;
  class checkValue;
  class macromolecule;
  class molecules;
  class salt;
  class grandcanonical;
    
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
   *
   *  <b>Displacement parameter optimization:</b>\n
   *  If dp_opt is set to true the displacement parameter, dp, will be optimized
   *  to give the highest possible mean square displacement (L^2). This is done by
   *  randomly generating a displacement parameter within the interval [dpmin:dpmax]
   *  and sample L^2 100 times. The resulting average will be added to a dp/L^2 table.
   *  When dp_N displacement parameters has been sampled, dp will be set to the value that gave the
   *  maximum L^2 and the optimization will automatically be turned off.
   */
  class markovmove {
    
  private:
    int dp_cnt;
    int dp_N;                                   //!< Number of displacement parameters to generate
    double dp_width;                            //!< Granularity of the distribution (same units as dp)
    double dp_min;                              //!< Minimum allowed displacement parameter
    double dp_max;                              //!< Maximum allowed displacement parameter
    xytable<double,average<double> > dp_dist;   //!< mean sq. displacement as a function of displacement parameter
    void init();                                //!< Initialization of variables (called from constructor)
    
  protected:
    slump slp;
    bool dp_opt;                                //!< True if displacement parameter should be optimized
    average<double> dpsqr;                      //!< Mean square displacement average
    double uold, unew, deltadp;
    double utot;                                //!< Sum of energy changes for all moves
    unsigned long long int cnt, naccept;
    string cite;                                //!< Reference to additional info. (article, url etc.)
    string prefix;                              //!< Inputfile keyword prefix
    container *con;
    ensemble *ens;
    bool run(float);                            //!< Probability
    virtual void getInput(inputfile &);         //!< Read parameters from inputfile
    virtual double newdp();                     //!< Generate new displacement parameter between [dp_min:dp_max]
    
  public:
    double optimaldp();                         //!< Retrieve optimal displacement parameter from DP/L^2 distribution
    string name;                                //!< Arbitrary name for the move
    enum keys {OK, ENERGY, HC};
    keys rc;                                    //!< Return code from move() functions
    double dp;                                  //!< Displacement parameter
    double du;                                  //!< Energy change of last move
    float runfraction;                          //!< Fractional chance that the move will be performed
    float accepted();                           //!< Return fraction of accepted moves
    void check(checkValue &);                   //!< Output test
    virtual string info();                      //!< Show info about group
    virtual double move();                      //!< Generic move function (should be called for every move)
    void adjust_dp(float=30, float=40);         //!< Adjust displacement parameter
    energybase *pot;
    markovmove(ensemble &, container &, energybase &);
    double zconstrain;                          //!< Constrains the displacement from getting closer than this value from a wall
  };
}//namespace
#endif
