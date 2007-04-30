#ifndef _markovmove_h
#define _markovmove_h

#include "container.h"
#include "potentials.h"
#include "ensemble.h"

/*! \brief Base class for MC moves
 *  \author Mikael Lund
 *  \todo Perhaps the P_pairpot could be made more elegant?
 *
 *  This class will keep track of the returned energy, if the
 *  move was successful of not and can provide statistics about the
 *  accepted and rejected moves.
 *
 *  Each derived class should provide a simple move() function that will
 *  perform a trial move, calculate the energy and either accept or reject
 *  the move. Unsuccessful moves are automatically undone.
 *
 *  The used pair-potential is identified using the
 *  type \verb T_pairpot \endverbatim that must be
 *  defined before processing the source code. For example,
 *
 *  \code
 *  #include "potentials.h"
 *  typedef pot_coulomb T_pairpot
 *  #include "markovmove.C"
 *  \endcode
 */
class markovmove {
  protected:
    double uold, unew;
    unsigned long long int cnt, naccept; 
    container *con;
    ensemble *ens;
    interaction<T_pairpot> *pot;
  public:
    enum keys {OK, ENERGY, HC};
    keys rc;
    double du,            //!< Energy change of last move
           utot;          //!< Sum of energy changes for all moves
    float accepted();     //!< Return fraction of accepted moves
    markovmove(ensemble &e, container &c, interaction<T_pairpot> &inter) {
      float du=utot=0;
      ens=&e;
      con=&c;
      pot=&inter;
    }
};

/*! \brief Move salt particles
 *  \author Mikael Lund
 */
class saltmove : public markovmove {
  public:
    saltmove( ensemble &, container&, interaction<T_pairpot>& );
    void move(group &, float, int=-1); //<! Move a salt particle
};

/*! \brief Symmetrically move two groups along z-axis
 *  \author Mikael Lund
 */
class dualzmove : public markovmove {
  public:
    float z;    //!< Distance between CM's of the groups
    void move(group &, group &);
};

/*! \brief Rotate group around mass-center.
*/
class rotate : public markovmove { 
  void move(group &);
};

/*! \brief Titrate all titrateable sites
 *  \author Mikael Lund
 */
class titrateall : public markovmove {
  //titrateall( space &, interaction<T_pairpot>&, container &, titrate &);
  void move(group &);
};

#endif
