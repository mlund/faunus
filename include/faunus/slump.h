#ifndef FAU_slump_h
#define FAU_slump_h
#include <cmath>
#include <cstdlib>
#include <ctime>

// Needed for sunCC
#include <stdlib.h>
#include <time.h>

#include "faunus/MersenneTwister.h"

namespace Faunus {
  class random {
    public:
      virtual double random_one()=0;              //!< Random number between [0:1[
      virtual void random_seed(unsigned int=0)=0; //!< Seed random generator (globally)
      bool runtest(float=0.5);                    //!< Probability bool
      double random_half();                       //!< Random number between [-0.5:0.5[
  };

  /*!
   * \brief Default C++ random number generator
   * \author mikaek lund
   * \date Lund, 2002
   */
  class randomDefault : public random {
    private:
      double rand_max_inv;
    public:
      randomDefault();
      void random_seed(unsigned int=0);
      double random_one();
  };

  /*!
   * \brief Mersenne Twister Random number functions
   * \author mikaek lund etc.
   * \date Prague, 2008
   */
  class randomTwister : public random {
    private:
      MTRand mt;
    public:
      void random_seed(unsigned int=0);
      double random_one();
  };

  /*! \typedef Speficy random number generator
   */
  typedef Faunus::randomDefault slump;
}
#endif

