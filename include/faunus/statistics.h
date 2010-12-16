#ifndef _statistics_h
#define _statistics_h
#include <vector>
#include <cmath>
namespace Faunus {
  /*!
   *  \brief Statistical functions for a single dataset
   *
   *  Can be used as the Faunus::average class but will store ALL
   *  added data in an internal vector so as to calculate
   *  rms, standard deviation etc.
   *
   *  Example:\n
   *  \code
   *  statistics<float> s;
   *  s+=2.
   *  s+=4.
   *  s.avg(); // --> 3.
   *  \endcode
   *
   *  \author Mikael Lund
   *  \date Prague, 2007
   *  \todo Specify weight for each sample. Via an internal vector?
   */
  template<class T>
    class statistics {
      private:
        std::vector<T> x;
      public:
        void operator+=(T); //!< Add a data point to set
        T avg();            //!< Return average of set
        T stddev();         //!< Return standard deviation
        T rms();            //!< Root-mean-square
    };

  template<class T>
    void statistics<T>::operator+=(T val) { x.push_back(val); }

  template<class T>
    T statistics<T>::avg()
    {
      T s=0;
      for (int i=0; i<x.size(); i++) s+=x[i];
      return s/x.size();
    }

  /*! \return \f$ s_{N-1} = \sqrt{\frac{1}{N-1} \sum \left ( x_i - \left <x\right > \right )^2}  \f$
  */
  template<class T>
    T statistics<T>::stddev()
    {
      int N=x.size();
      if (N<2) return 0;
      T xm=avg();
      T sum=0;
      for (int i=0; i<N; i++)
        sum+=pow( x[i]-xm, 2 );
      return sqrt( sum / (N-1) );
    }

  /*! \return \f$ rms = \sqrt{\left <x^2 \right> }\f$
  */
  template<class T>
    T statistics<T>::rms()
    {
      T sum=0;
      for (int i=0; i<x.size(); i++)
        sum+=x[i]*x[i];
      return sqrt(sum/x.size());
    }
}
#endif
