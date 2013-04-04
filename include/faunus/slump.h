#ifndef FAU_slump_h
#define FAU_slump_h

#include <string>
#include <random>
#include <cassert>

namespace Faunus {

  /*@
   * @brief Base class for random number generation
   *
   * Derived classes need only provide `_randone()` and `seed()` function.
   */
  class RandomBase {
    private:
      virtual double _randone()=0;  //!< Random number in range `[0:1[`
    public:
      std::string name;
      virtual ~RandomBase();
      virtual void seed(int=0)=0;   //!< Seed random generator
      double randHalf();            //!< Random number in range `[-0.5:0.5[`
      unsigned int rand();          //!< Random number in range `[0:max unsigned int[`
      /*!
       * \brief Random number in range [0:1[
       */
      inline double operator()() {
        double x=_randone();
        assert(x>=0 && x<1 && "Random number out of range!");
        return x;
      }
  };

 /**
  * @brief Ran2 Random Number Gererator
  * @author Bjorn Persson
  * @date Lund, 2008
  * @note A class for ran2 from 'Numerical Recipies'.
  * @warning Not thread safe!
  */
  class RandomRan2: public RandomBase {
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
      double _randone();
    public:
      RandomRan2();
      void seed(int=-7);
  };

  /**
   * @brief Mersenne Twister Random number functions (C++11)
   *
   * This uses the C++11 Mersenne Twister for random number generations and while
   * quite slow ran2, MT generally provides better randomness.
   *
   * @date Lund, 2010
   */
  template<typename T=double, typename Tengine=std::mt19937> 
    class RandomTwister : public RandomBase {
      private:
        Tengine eng; //pseudo random number engine
        std::uniform_real_distribution<T> dist;
        T _randone() {
          T x;
#pragma omp critical
          x=dist(eng);
          return x;
        }
      public:
        RandomTwister() : dist(0,1) {
          name="Mersenne Twister";
          //std::random_device rd;
          //eng.seed( rd() );
        }
        void seed(int s) {
#pragma omp critical
          eng.seed(s);
        }
    };

#if defined(MERSENNETWISTER)
  typedef Faunus::RandomTwister<double,std::mt19937> slump;
#else
  typedef Faunus::RandomRan2 slump;
#endif
  extern slump slp_global;
}
#endif
