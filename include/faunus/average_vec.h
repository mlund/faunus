#ifndef FAU_average_vec_h
#define FAU_average_vec_h

#include <iomanip>
#include <numeric>

namespace Faunus {
  /*
   * Average class (BIOSIM VERSION)
   *
   * Convenient handling of average
   * values.
   *
   * M. Lund, 2002-2007
   *
   */

  //! Class to collect average values.
  /*!
   * This will automatically count the number of
   * values you add and can even return the standard
   * deviation (assuming the values are independent).
   * Example:\n
   *   average x;\n
   *   x+=2.0;\n
   *   x+=6.0;\n
   *   cout << x ;   // --> 4.0\n
   */

  template<class T>
    class average_vec {
      public:
      std::vector<T> v;                    ///< Vector to store sets of length average::maxcnt
      unsigned int maxcnt;            ///< Number of values in each set.
      unsigned long long int cnt;     ///< Number of values in current set.
      average_vec(int=10000);             ///< Constructor. Specify average::maxcnt
      T sum;                          ///< Sum of current set.
      T operator+(average_vec &);      
      T operator-(average_vec &);      
      void operator+=(T);             ///< Add value to current set. Same is average::add
      void add(T);                    ///< Add value to current set.
      void info();
      void reset();                   ///< Clear all data
      T avg();                        ///< Return average of all sets in average::v
      T stdev();                      ///< If more than one set, return standard deviation.
    };

  template<class T>
    average_vec<T>::average_vec(int max) {
      maxcnt=max;
      sum=cnt=0;
      v.clear();
    }

  template<class T>
    void average_vec<T>::reset() {
      sum=cnt=0;
      v.clear();
    }

  template<class T>
    T average_vec<T>::operator+(average_vec &a) { return avg() + a.avg(); };

  template<class T>
    T average_vec<T>::operator-(average_vec &a) { return avg() - a.avg(); };

  template<class T>
    void average_vec<T>::operator+=(T x) { add(x); };

  template<class T>
    void average_vec<T>::add(T x) {
      sum+=x;
      cnt++;
      if (cnt==maxcnt) {
        v.push_back( sum/cnt );
        sum=0;
        cnt=0;
      }
    }

  template<class T>
    T average_vec<T>::avg() {
      return (v.size()==0) ? sum/cnt :
        std::accumulate(v.begin(),v.end(),0.)/v.size();
    }

  template<class T>
    T average_vec<T>::stdev() {
      T vav=avg(), sum=0;
      for (unsigned int i=0; i<v.size(); i++)
        sum+=pow(v[i]-vav, 2);
      return sqrt( sum / ( v.size() * (v.size()-1) ) );  
    }

  template<class T>
    void average_vec<T>::info() {
      std::cout << "# " << v.size() << " sets each containing " << maxcnt << " points." << std::endl;
    }
}
#endif
