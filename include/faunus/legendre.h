#ifndef FAU_LEGENDRE_H
#define FAU_LEGENDRE_H

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

