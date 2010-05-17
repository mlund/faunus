#ifndef FAU_slump_h
#define FAU_slump_h
#include <cmath>
#include <iostream>
#include <sstream>
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
    protected:
      std::string name;
    public:
      virtual double random_one()=0;              //!< Random number between [0:1[
      virtual void random_seed(int=0)=0;          //!< Seed random generator (globally)
      bool runtest(float=0.5);                    //!< Probability bool
      double random_half();                       //!< Random number between [-0.5:0.5[
      std::string info();                         //!< Print information string
  };

  /*!
   * \brief Default C++ random number generator
   * \author Mikael Lund
   * \date Lund, 2002
   * \warning random_one sometimes returns 1 (one)!!
   */
  class randomDefault : public random {
    private:
      double rand_max_inv;
    public:
      randomDefault();
      void random_seed(int=0);
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
      randomTwister();
      void random_seed(int=0);
      double random_one();
  };


 /*!
  * \brief Ran2 Random Number Gererator
  * \author Bjorn Persson
  * \date Lund, 2008
  * \note A class for ran2 from 'Numerical Recipies'.
  * \warning Not thread safe!
  */
  class ran2: public random {
    private:
      const static int IM1=2147483563, IM2=2147483399;
      const static int IA1=40014, IA2=40692, IQ1=53668, IQ2=52774;
      const static int IR1=12211, IR2=3791, NTAB=32, IMM1=IM1-1;
      const static int NDIV=1+IMM1/NTAB;
      const static double EPS;//=3.0e-16;
      double AM,RNMX;
      int idum, idum2, iy;
      int iv[32]; 
      int j,k;
      double temp;
    public:
      ran2();
      double random_one();
      void   random_seed(int=-7);
  };

  //typedef Faunus::randomDefault slump; // Generator selection!
  typedef Faunus::ran2 slump; // Generator selection!
  
  extern slump slp;

}
#endif
