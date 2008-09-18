/*! \brief Integration using Simpson's method
 *  \author mikaek lund 
 *  \param a Lower boundary
 *  \param b Upper boundary
 *  \param n Number of steps
 *
 *  \example simpson-example.C
 */
template<double f(double)>
double simpson(double a,double b, int n) {
  double step,r,res;
  step=(b-a) / 2. / n;
  r=f(a);
  res=(r + f(b))/2.;
  for (int i=1; i<2*n; i++) {
    r=f( a + i*step );
    if ( (i%2)!=0 )
      res += r+r;
    else res += r; 
  }
  res *= step*2./3.;
  return res;
}
