#pragma once

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
          return std::sqrt( (sqsum + cnt*avg()*avg() - 2*sum*avg()) / static_cast<double>(cnt-1) );
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

#ifdef DOCTEST_LIBRARY_INCLUDED
    TEST_CASE("[Faunus] Average") {
        Average<double> a;
        a+=1.0;
        a+=2.0;
        CHECK( a.cnt == 2 );
        CHECK( a.sum == 3 );
        CHECK( a.sqsum == 5 );
        CHECK( a.avg() == 1.5 );
        CHECK( a == 1.5 ); // implicit conversion to double

        auto b = a; // copy
        CHECK( !b.empty() ); // check not empty()
        CHECK( a==b );
        b.clear(); // reset all data
        CHECK( b.empty() ); // check empty()
        b+=2.0;
        b+=3.0;
        CHECK( a<b ); // a.avg() < b.avg()
        CHECK( (a+b).avg() == doctest::Approx(2) );

        b = 1.0; // assign from double
        CHECK( b.size()==1 );
    }
#endif

}//namespace
