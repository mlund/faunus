#undef NDEBUG
#include <faunus/faunus.h>

using namespace Faunus;

bool eq(double a, double b, double tol=1e-6) { return (std::abs(a-b)<tol) ? true : false; }

int main() {

  assert( pc::infty == -std::log(0) && "Problem with infinity" );

  assert( eq( exp_cawley(1),std::exp(1),1e-1 ) && "Problem with approximate exp() function" );

  // check geometries
  Geometry::Sphere geoSphere(1000);
  Geometry::Cylinder geoCylinder(1000,1000);
  Point a( 0,0,sqrt(16)), b(0,sqrt(64),0);
  double x = geoSphere.sqdist(a,b);
  double y = geoCylinder.sqdist(a,b);
  assert( x==16+64 && y==16+64 );
  assert( eq(x,y) && "No good distance calculation");

  // check table of averages
  typedef Analysis::Table2D<float,Average<float> > Ttable;
  Ttable table(0.1, Ttable::XYDATA);
  table(2.1)+=1;
  table(2.1)+=3;
  assert( eq( table(2.1), 2 ) && "Bad average or 2D table");
}
