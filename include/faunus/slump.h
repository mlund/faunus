#ifndef FAU_slump_h
#define FAU_slump_h
#include <cmath>
#include <cstdlib>  // Including this SUCKS! Be aware of abs()!
#include <ctime>

// Needed for sunCC ("__SUNPRO_CC" or "__sun" ?)
#ifdef __SUNPRO_CC
#include <stdlib.h>
#include <time.h>
#endif

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
   * \author Mikael Lund
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
   * \author Mikael Lund
   * \date Prague, 2008
   * \note Not thread safe
   * \warning Weird behavior with rotational moves.
   * \todo Fix weird behavior
   *
   * To enable this random number generator -DTWISTER must be
   * set at compile time.
   */
  class randomTwister : public random {
    private:
      MTRand mt;
    public:
      void random_seed(unsigned int=0);
      double random_one();
  };

  typedef Faunus::randomDefault slump; // Generator selection!

}
#endif
