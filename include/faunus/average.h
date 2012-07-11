#ifndef FAU_AVERAGE_H
#define FAU_AVERAGE_H
//%template(average_dbl) average<double>;

#ifndef SWIG
#include "faunus/common.h"
#endif

namespace Faunus {
  /*!
   * \brief Class to collect average values
   * \author Mikael Lund
   * \date 2002-2010
   *
   * \code
   * average<double> x,y;
   * x+=2.0;
   * x+=6.0;
   * y+=1.0;
   * double f=x+y;
   * cout << x << " " << y << " " << f; // --> 4.0 1.0 3.0 
   * \endcode
   */
  template<class T> class Average {
  protected:
    //   vector v;
  public:
    T sqsum;                                      ///< Square sum
    Average();                                 
    virtual ~Average() {} 
    T sum;                                        ///< Sum of all values
    unsigned long long int cnt;                   ///< Number of values
    T avg() const;                                ///< Return average
    T rms();                                      ///< Root-mean-square
    T stdev();                                    ///< Standard deviation
    void add(T);                                  ///< Add value to current set.
    void reset();                                 ///< Clear all data
    Average & operator=(T);                       ///< Assign value to current set. 
    Average & operator+=(T);                      ///< Add value to current set. 
    Average & operator*=(T);                      ///< Scale current set
    Average & operator/=(T);                      ///< Scale current set
    operator T() const;                           ///< Static cast operator
    T operator*(T) const;                         ///< Evaluates average times T
    T operator+(T) const;                         ///< Evaluates average plus T
    const Average operator+(const Average&) const;///< Merge two averages (with correct weights)
    bool operator==(const Average &) const;       ///< Comparison operator
    bool operator<(const Average &) const;
  };
  
  template<class T> Average<T>::Average() { reset(); }
  
  template<class T> T Average<T>::avg() const {
    assert(cnt>0 && "Average counter is zero");
    return double(sum)/cnt;
  }

  template<class T> void Average<T>::reset() {
    sum=sqsum=0;
    cnt=0;
  }

  template<class T> Average<T>& Average<T>::operator=(T x) {
    reset();
    add(x);
    return *this;
  }

  template<class T> Average<T>::operator T() const { return avg(); }

  template<class T> T Average<T>::operator*(T x) const { return avg() * x; }

  template<class T> T Average<T>::operator+(T x) const { return avg() + x; }
  
  template<class T> const Average<T> Average<T>::operator+(const Average &a) const {
    Average<T> r = *this;
    r.cnt += a.cnt;
    r.sum += a.sum;
    r.sqsum+=a.sqsum;
    return r;
  }

  /*!
   * \note sqsum is not scaled.
   */
  template<class T> Average<T> & Average<T>::operator*=(T x) {
    sum*=x;
    return *this;
  }

  /*!
   * \note sqsum is not scaled.
   */
  template<class T> Average<T> & Average<T>::operator/=(T x) {
    sum/=x;
    return *this;
  }
 
  template<class T> Average<T> & Average<T>::operator+=(T x) {
    add(x);
    return *this;
  }
  
  template<class T> void Average<T>::add(T x) {
    assert( cnt+1>0 && "Counter has rotated to zero!");
    if (cnt+1<=0)
      std::cerr << "Numeric limits reached in average class!\n";
    else {
      sum+=x;
      sqsum+=x*x;
      cnt++;
    }
  }

  template<class T> bool Average<T>::operator==(const Average &a) const {
    return (*this==a);
  }

  template<class T> bool Average<T>::operator < (const Average &a) const {
    //if (cnt==0)
    //  return true;
    return avg() < a.avg();
  }

  template<class T> T Average<T>::rms() {
    return sqrt(sqsum/cnt);
  }

  template<class T> T Average<T>::stdev() {
    return sqrt( sqsum/cnt - pow(sum/cnt,2) );
  }

  template<class T> class AverageExt {
    private:
      long unsigned int nmax, cnt;
      std::list<T> v;
    public:
      AverageExt(unsigned int maxNumberOfSamples=1e9) : v(1) {
        nmax = maxNumberOfSamples;
        v.back()=0;
        cnt=0;
      }
      AverageExt& operator+=(T x) {
        cnt++;
        v.back()+=x;
        if (cnt==nmax) {
          v.back() /= nmax;
          v.push_back(0);
          cnt=0;
        }
        return *this;
      }

      T avg() const {
        if (v.size()>1)
          return std::accumulate(v.begin(),v.end(),0.) / v.size();
        return (cnt==0) ? 0 : v.back()/cnt;
      }

      T stddev() const {
        if (v.size()>1) {
          T sum=0, vav=avg();
          for (auto i : v)
            sum+=pow( i-vav, 2 );
          return sqrt( sum / (v.size()-1) );
        }
        return 0;
      }

      friend std::ostream &operator<<(std::ostream &o, const AverageExt<T> &a) {
        double s=a.stddev();
        o << a.avg();
        if (s>1e-15)
          o << " " << s;
        return o;
      } 
  };

  /*!
   * \brief Class to keep track of block correlations
   * \author Mikael Lund
   * \date October 2009
   *
   * The sampling is performed in blocks of length n specified in
   * the constructor.
   *
   * \f$ c_i = \frac{ \langle x_0x_i\rangle_{i<n} - \langle x\rangle^2 }
   * { \langle x^2\rangle - \langle x\rangle^2  } \f$
   *
   * Example:
   * \code
   * BlockCorrelation ci(50); // energy correlation
   * ci += uinit + du;        // place in MC loop
   * ...
   * for (size_t i=0; i<ci.size(); i++)
   *   cout << i << " " << ci[i] << end;
   * \endcode
   * ci will eventually fall off from one (full correlation) to
   * zero (uncorrelated).
   */
  template<class T=double>
    class BlockCorrelation {
      private:
        unsigned int n,            //!< Length of each correlation measurement
                     cnt;          //!< Internal counter for each correlation set
        vector< Average<T> > xixj; //!< Average correlation product, <xixj>
        Average<T> xmean;          //!< Average values, <x>
        T xi;                      //!< Reference value (i=0) for each correlation set
      public:
        BlockCorrelation(unsigned int=50);
        BlockCorrelation & operator+=(T);    //!< Sample value
        T operator[] (unsigned int);    //!< Get correlation at i
        unsigned int size();            //!< Get block length
        bool dump(string filename); //!< Dump to disk
    };

  /*!
   * \param len Sample length
   */
  template<class T>
    BlockCorrelation<T>::BlockCorrelation(unsigned int len) : cnt(0), n(len) {
      xixj.resize(len);
    }

  template<class T>
    BlockCorrelation<T> & BlockCorrelation<T>::operator+=(T xj) {
      xmean+=xj;            // update total average
      if (cnt==n)           // end of block? reset counter.
        cnt=0;
      if (cnt==0)           // block start? save reference value.
        xi=xj;
      xixj.at(cnt) += xi*xj;// update average
      cnt++;
      return *this;
    }

  template<class T>
    T BlockCorrelation<T>::operator[] (unsigned int i) {
      T xm=xmean.avg();
      T x2m=xmean.sqsum/xmean.cnt;
      return ( xixj.at(i).avg() - xm*xm ) / ( x2m - xm*xm ); 
    }

  template<class T>
    unsigned int BlockCorrelation<T>::size() {
      return xixj.size();
    }

  //! Dump to disk
  template<class T>
    bool BlockCorrelation<T>::dump(string filename) {
      std::ofstream f(filename.c_str());
      if (f) {
        f.precision(6);
        for (int i=0; i<size(); i++)
          f << operator[](i) << endl; 
        return true;
      }
      return false;
    }

}//namespace
#endif
