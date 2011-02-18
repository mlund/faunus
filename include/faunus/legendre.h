#ifndef FAU_LEGENDRE_H
#define FAU_LEGENDRE_H

namespace Faunus {
  namespace Math {
#ifdef FAU_FAST_SQRT
    /*!
     * \brief Log Base 2 Approximation for the square root
     * \note http://ilab.usc.edu/wiki/index.php/Fast_Square_Root
     * 
     * Roughly 5 times faster than sqrt() - four per percent error.
     */
    inline double fastSqrt_2(const double x) {
      union {
        long i;
        double x;
      } u;
      u.x = x;
      u.i = (((long)1)<<61) + (u.i >> 1) - (((long)1)<<51); 
      return u.x;
    }

    /*!
     * \brief Log Base 2 Approximation With One Extra Babylonian Steps for sqrt()
     * \note http://ilab.usc.edu/wiki/index.php/Fast_Square_Root
     *
     * Approximately two times faster than sqrt()
     */
    inline double fastSqrt_Bab(const double x) {
      union {
        long i;
        double x;
      } u;
      u.x = x;
      u.i = (((long)1)<<61) + (u.i >> 1) - (((long)1)<<51); 
      u.x = 0.5F * (u.x + x/u.x);
      return u.x;
    }

    #define SQRT_MAGIC_D 0x5fe6ec85e7de30da
    /*!
     * \brief Quake 3 Walsh Method for the inverse square root
     */
    inline double invSqrt_Q3(const double x) {
      const double xhalf = 0.5F*x;
      union {
        double x;
        long i;
      } u;
      u.x = x;
      u.i = SQRT_MAGIC_D - (u.i >> 1);  // gives initial guess y0
      return u.x*(1.5F - xhalf*u.x*u.x);// Newton step, repeating increases accuracy
    }
    inline double fastSqrt_Q3(const double x) {
      return x * invSqrt_Q3(x);
    }
#endif

  } //math namespace
} //faunus namespace

namespace Faunus {
  /*!
   * \brief Evaluate n'th degree Legendre polynomium
   * \author Mikael Lund
   * \date Canberra 2005-2006
   *
   * Example\n
   * \code
   * legendre l(10);
   * l.eval(1.3);
   * cout << l.p[3]
   * \endcode
   */
  class legendre {
    private:
      unsigned short n;           //!< Legendre order
    public:
      std::vector<double> p;      //!< Legendre terms stored here
      //! Evaluate polynomium at x
      inline void eval(double x) {
        if (n > 0) {
          p[1] = x;
          double di;
          for( unsigned short i=1; i<n; i++ ) {
            di=static_cast<double>(i);
            p[i+1] = (2*di + 1) / (di+1) * x * p[i] - di/(di+1)*p[i-1];
          }
        }
      }

      /*!
       * \param order Order of the polynomium
       */
      legendre(unsigned short order=0) { resize(order); }

      void resize(unsigned short order) {
        n=order;
        p.resize(n+1);
        p[0]=1.;
      }
  };
}//namespace
#endif

