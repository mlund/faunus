#include<vector>
#include<iostream>

using namespace std;

/*!
 * \brief Class to handle a 2D table.
 * \author Mikael Lund
 */
class table{
  private:
    vector<double> y;
    double min,max,steps;
    unsigned short x2i(double x)
    { return static_cast<unsigned short>(x/steps+.5); }
  public:
    table( double beg, double end, double resolution) {
      min=beg;
      max=end;
      steps=resolution;
      y.resize( x2i(max)+1 );
    }
    void set(double x, double val) {
      if (x>=min && x<=max)
        y[ x2i(val) ]=val;
    }
    double get(double x) {
      if (x>=min && x<=max)
        return y[ x2i(x) ];
    }
};

/*!
 * \brief Evaluate n'th degree Legendre polynomium at x
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
    void resize(unsigned short);//!< Re-set order
 public:
    vector<double> p;           //!< Legendre terms stored here
    legendre(unsigned short); 

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
 
