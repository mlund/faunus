#undef NDEBUG
#include <faunus/faunus.h>

using namespace Faunus;

int main() {

  assert( pc::infty == -std::log(0) && "Problem with infinity" );

  assert( std::abs( exp_cawley(1)-std::exp(1) )<1e-1 && "Problem with approximate exp() function" );

  // check geometries
  Geometry::Sphere geoSphere(1000);
  Geometry::Cylinder geoCylinder(1000,1000);
  Point a( 0,0,sqrt(16)), b(0,sqrt(64),0);
  double x = geoSphere.sqdist(a,b);
  double y = geoCylinder.sqdist(a,b);
  assert( x==16+64 && y==16+64 );
  assert( abs(x-y)<1e-6 && "No good distance calculation");
}
