#ifndef FAU_slump_h
#define FAU_slump_h

#ifdef SWIG
%module pyfaunus
%{
#include "faunus/slump.h"
%}
#else
  #include <cmath>
  #include <iostream>
  #include <sstream>
  #include <cstdlib>  // Including this SUCKS! Be aware of abs()!
  #include <ctime>
  #include <random>

  // Needed for sunCC ("__SUNPRO_CC" or "__sun" ?)
  #ifdef __SUNPRO_CC
    #include <stdlib.h>
    #include <time.h>
  #endif
#endif

namespace Faunus {
  class random {
    protected:
      std::string name;
    public:
      virtual ~random() {}
      virtual double random_one()=0;              //!< Random number between [0:1[
      virtual void random_seed(int=0)=0;          //!< Seed random generator (globally)
      bool runtest(float=0.5);                    //!< Probability bool
      double random_half();                       //!< Random number between [-0.5:0.5[
      std::string info();                         //!< Print information string
      virtual unsigned int rand()=0;              //!< Random number between 0 and rand_max
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
      unsigned int rand();
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
      static const int IM1=2147483563, IM2=2147483399;
      static const int IA1=40014, IA2=40692, IQ1=53668, IQ2=52774;
      static const int IR1=12211, IR2=3791, NTAB=32, IMM1=IM1-1;
      static const int NDIV=1+IMM1/NTAB;
      static const double EPS;//=3.0e-16;
      double AM,RNMX;
      int idum, idum2, iy;
      int iv[32]; 
      int j,k;
      double temp;
    public:
      ran2();
      double random_one();
      void   random_seed(int=-7);
      unsigned int rand();
  };

  /*!
   * \brief Mersenne Twister Random number functions (C++ TR1)
   * \author Mikael Lund
   * \date Lund, 2010
   */
  class randomTwister : public random {
    private:
      double maxinv;
      std::mt19937 eng;
      std::uniform_real_distribution<double> dist;
    public:
      randomTwister();
      double random_one();
      void random_seed(int=0);
      unsigned int rand();
  };

#if defined(MERSENNETWISTER)
  typedef Faunus::randomTwister slump;
#else
  typedef Faunus::ran2 slump;
#endif
  extern slump slp_global;
}
#endif
