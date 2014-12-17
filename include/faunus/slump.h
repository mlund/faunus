#ifndef FAU_slump_h
#define FAU_slump_h

#include <string>
#include <random>
#include <cassert>

/// @brief Namespace for Faunus
namespace Faunus {

  /**
   * @brief Mersenne Twister Random number generator
   * @date Lund, 2010
   */
  template<typename T=double, typename Tengine=std::mt19937> 
    class RandomTwister {

      private:

        Tengine eng; //random number engine

        std::uniform_real_distribution<T> dist;

        T _randone() {
          T x;
#pragma omp critical
          x=dist(eng);
          return x;
        }

      public:

        /// @brief Constructor -- default deterministic seed is used
        RandomTwister() : dist(0,1) {}

        /// @brief Integer in range [min:max]
        int range(int min, int max) {
          std::uniform_int_distribution<int> d(min,max);
          return d(eng);
        }

        /**
         * @brief Seed random number engine
         *
         * The default value is obtained from `std::random_device`
         * which will give a non-deterministic run.
         * An 32-bit integer can be given as well for a
         * deterministic seed. Note that `seed()`
         * is *not* called upon construction.
         */
        template<typename Tint>
          void seed(Tint s=std::random_device()) {
#pragma omp critical
            eng.seed(s);
          }

        /// @brief Random number in range `[-0.5,0.5)`
        T randHalf() { return _randone() - 0.5; }

        /// @brief Random number in range [0,1)
        T operator()() { return _randone(); }

        /**
         * @brief Random number in range `[0,max unsigned int)`
         * @todo To be removed
         */
        unsigned int rand() {
          return range(0, std::numeric_limits<unsigned int>::max()-1);
        }
    };

  typedef Faunus::RandomTwister<double,std::mt19937> slump;

  extern slump slp_global;

  /**
   * @brief Get random element from iterable container (vector, map, group etc.)
   * @return Iterator to random element
   */
  template<class Titer>
    Titer randomElement(const Titer &beg, const Titer &end) {
      auto i=beg;
      std::advance(i, slp_global.range( 0,std::distance(beg, end)-1 ) );
      return i;
    }

} // namespace

#endif
