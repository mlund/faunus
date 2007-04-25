#ifndef _legendre_h
#define _legendre_h
#include<vector>

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
    legendre(unsigned short=0); 
    void resize(unsigned short);//!< Re-set order
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
};

/*!
 * \param order Order of the polynomium
 */
legendre::legendre(unsigned short order) { resize(order); }
void legendre::resize(unsigned short order) {
  n=order;
  p.resize(n+1);
  p[0]=1.;
}

#endif

