#pragma once
#include <doctest/doctest.h>
#include <limits>
#include <ostream>
#include <istream>
#include <cmath>

namespace Faunus
{

  /** @brief Class to collect averages */
  template<class T=double> struct Average
  {
      T sqsum=0;                    //!< Square sum
      unsigned long long int cnt=0; //!< Number of values
      T sum=0;                      //!< Sum

      template <class Archive> void serialize(Archive &archive) { archive(sum, sqsum, cnt); }

      double avg() const {
          return sum / static_cast<double>(cnt);
      } //!< Average

      double rms() const {
          return empty() ? 0 : std::sqrt(sqsum / cnt);
      } //!< Root-mean-square

      double stdev() const {
          return std::sqrt((sqsum + static_cast<T>(cnt) * avg() * avg() - 2 * sum * avg()) /
                           static_cast<double>(cnt - 1));
      } //!< Standard deviation

      void add(T x) {
          if ( cnt+1==0 || sqsum + x*x > std::numeric_limits<T>::max() )
              throw std::runtime_error("Numeric overflow");
          cnt++;
          sum += x;
          sqsum += x*x;
      } //!< Add value to current set

      void clear() {
          sum = sqsum = 0;
          cnt = 0;
      } //!< Clear all data

      bool empty() const {
          return cnt==0;
      } //!< True if empty

      auto size() const { return cnt; } //!< Number of samples

      Average& operator=(T x) {
          clear();
          add(x);
          return *this;
      } //!< Assign value to current set

      Average& operator+=(T x) {
          add(x);
          return *this;
      } //!< Add value to current set

      operator double() const {
          return avg();
      } //!< Static cast operator

      const Average operator+( const Average &a ) const {
          Average<T> r = *this;
          r.cnt += a.cnt;
          r.sum += a.sum;
          r.sqsum += a.sqsum;
          return r;
      } //!< Merge two averages (with correct weights)

      bool operator<( const Average &a ) const { return avg() < a.avg(); } //!< Compare operator

      friend std::ostream &operator<<( std::ostream &o, const Average<T> &a ) {
          o << a.cnt << " " << a.sum << " " << a.sqsum;
          return o;
      } // serialize to stream

      Average<T> &operator<<( std::istream &in ) {
          in >> cnt >> sum >> sqsum;
          return *this;
      } // de-serialize from stream
  };

#if __cplusplus > 201703L
  /** Requirements for types used with `AverageObj` */
  template <class T> concept Averageable = requires(T a) {
      a * 1.0; // please implement `T operator*(double) const`
      a *a;    // please implement `T operator*(const T&) const`
      a += a;  // please implement `T& operator+=(const T&)`
  };
#endif
      /**
       * @brief Simple class to average data contained in objects
       * @tparam T Type to average
       * @tparam int_t Unsigned interger type
       *
       * It is required that `T` has the following operator overloads:
       * - `T operator*(double) const`
       * - `T operator*(const T&) const`
       * - `T& operator+=(const &T)`
       */
#if __cplusplus > 201703L
  template <Averageable T, typename int_t = unsigned long int> class AverageObj {
#else
  template <typename T, typename int_t = unsigned long int> class AverageObj {
#endif
    protected:
      int_t number_of_samples = 0;
      T sum; // make sure constructors zero this!
    public:
      AverageObj() : sum(T()){}; //!< Construct from empty object

      AverageObj(const T &value) : number_of_samples(1), sum(value){};

      //! Add to average
      AverageObj &operator+=(const T &value) {
          if (number_of_samples == std::numeric_limits<int_t>::max()) {
              throw std::overflow_error("maximum samples reached");
          } else {
              sum += value;
              ++number_of_samples;
          }
          return *this;
      }
      //! Calculate average
      T avg() const {
          if (number_of_samples > 0) {
              return sum * (1.0 / static_cast<double>(number_of_samples));
          } else {
              return sum;
          }
      }

      //! Convert to T
      operator T() const { return avg(); }

      //! Compare operator
      bool operator<(const AverageObj &other) const { return avg() < other.avg(); }

      //! True if empty
      bool empty() const { return number_of_samples == 0; }

      //!< Clear all data
      void clear() { *this = AverageObj<T>(); }

      //! Number of samples
      int_t size() const { return number_of_samples; }
  };

  /*
   * Experimental extension of AverageObj that includes stdev() and rms(). This
   * requires a pow(T, double) function for squaring and taking the square-root.
   * How could this best be implemented and maintain compatibility when T=double?
   */
  template <typename T, typename int_t = unsigned long int> class AverageObjStdev : public AverageObj<T> {
    private:
      T sum_squared; // make sure constructors zero this!
      using AverageObj<T>::avg;
      using AverageObj<T>::sum;
      using AverageObj<T>::number_of_samples;

    public:
      AverageObjStdev() : sum_squared(T()){}; //!< Construct from empty object

      AverageObjStdev(const T &value) : AverageObj<T>(value), sum_squared(value * value){};

      //! Add to average
      AverageObjStdev &operator+=(const T &value) {
          AverageObj<T>::operator+=()(value);
          sum_squared += std::pow(value, 2);
          return *this;
      }

      //! Root-mean-square
      T rms() const {
          if (number_of_samples == 0) {
              return T();
          } else {
              return std::pow(sum_squared * (1.0 / static_cast<double>(number_of_samples)), 0.5);
          }
      }

      //! Standard deviation
      T stdev() const {
          if (number_of_samples == 0) {
              return T();
          } else {
              const auto N = static_cast<double>(number_of_samples);
              const auto mean = avg();
              return std::pow((sum_squared + N * mean * mean - 2.0 * sum * mean) * (1.0 / (N - 1.0)), 0.5);
          }
      }

      //!< Clear all data
      void clear() { *this = AverageObjStdev<T>(); }
  };

  } // namespace Faunus
