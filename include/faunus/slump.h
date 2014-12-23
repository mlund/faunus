#ifndef FAU_slump_h
#define FAU_slump_h

#include <random>

/// @brief Namespace for Faunus
namespace Faunus {

  /**
   * @brief Mersenne Twister Random number generator for uniform distribution
   * @date Lund, 2010
   *
   * By default this generator is seeded with a default seed,
   * resulting in a deterministic sequence. For attempting
   * a non-deternimistic sequence, call `seed()` with no
   * arguments to activate a hardware induced seed using
   * `std::random_device`.
   */
  template<typename T=double, typename Tengine=std::mt19937>
    class RandomTwister {

      private:
        Tengine eng; //random number engine
        std::uniform_real_distribution<T> dist;

      public:

        /** @brief Constructor -- default deterministic seed */
        RandomTwister() : dist(0,1) {}

        /** @brief Integer in uniform range [min:max] */
        int range(int min, int max) {
          std::uniform_int_distribution<int> d(min,max);
          return d(eng);
        }

        /** @brief Seed random number engine (s>0: deterministic, s=0 non-deterministic) */
        void seed(unsigned int s=0) {
          if (s==0) {
#pragma omp critical
            eng = Tengine( std::random_device()() );
          } else {
#pragma omp critical
            eng.seed(s);
          }
        }

        /** @brief Random number in uniform range [0,1) */
        T operator()() {
          T x;
#pragma omp critical
          x=dist(eng);
          return x;
        }

        /** @brief Random number in uniform range `[-0.5,0.5)` */
        T randHalf() { return operator()() - 0.5; }
    };

  extern RandomTwister<> slp_global;

  /**
   * @brief Get random element from iterable container (vector, map, group etc.)
   * @return Iterator to random element
   * @note Constant-time complexity for random access iterators; otherwise linear.
   */
  template<class Titer>
    Titer randomElement(const Titer &beg, const Titer &end) {
      auto i=beg;
      std::advance(i, slp_global.range( 0, std::distance(beg, end)-1 ) );
      return i;
    }

} // namespace

#endif
