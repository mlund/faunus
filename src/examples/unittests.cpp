#undef NDEBUG
#include <faunus/faunus.h>

using namespace Faunus;

bool eq(double a, double b, double tol=1e-6) { return (std::abs(a-b)<tol) ? true : false; }

int main() {

  // check infinity
  assert( pc::infty == -std::log(0) && "Problem with infinity" );

  // check approximate exp()
  assert( eq( exp_cawley(1),std::exp(1),1e-1 ) && "Problem with approximate exp() function" );

  // check approximate 1/sqrt()
  assert( eq( invsqrtQuake(20.), 1/std::sqrt(20.), 1e-1) && "Problem w. inverse Quake sqrt.");

  // check geometries
  Geometry::Sphere geoSph(1000);
  Geometry::Cylinder geoCyl(1000,1000);
  Point a( 0,0,sqrt(16)), b(0,sqrt(64),0);
  double x = geoSph.sqdist(a,b);
  double y = geoCyl.sqdist(a,b);
  assert( x==16+64 && y==16+64 );
  assert( eq(x,y) && "No good distance calculation");

  // check random numbers
  {
    int min=10, max=0, N=1e7;
    double x=0;
    for (int i=0; i<N; i++) {
      int j = slp_global.rand() % 10;
      if (j<min) min=j;
      if (j>max) max=j;
      x+=j;
    }
    assert( min==0 && max==9 );
    assert( std::abs(x/N-4.5) < 1e-2 ); // average should be 4.5
  }

  // check vector rotation
  {
    Geometry::VectorRotate vrot;
    a.clear();
    a.x()=1.;
    vrot.setAxis( geoCyl, Point(0,0,0), Point(0,1,0), pc::pi/2); // rotate around y-axis
    a = vrot(a); // rot. 90 deg.
    assert( eq(a.x(),0,1e-8) && "Vector rotation failed");
    a = vrot(a); // rot. 90 deg.
    assert( eq(a.x(),-1,1e-8) && "Vector rotation failed");
  }

  {
    Geometry::QuaternionRotate qrot;
    a.clear();
    a.x()=1.;
    qrot.setAxis( geoCyl, Point(0,0,0), Point(0,1,0), pc::pi/2); // rotate around y-axis
    a = qrot(a); // rot. 90 deg.
    assert( eq(a.x(),0,1e-8) && "Vector rotation failed");
    a = qrot(a); // rot. 90 deg.
    assert( eq(a.x(),-1,1e-8) && "Vector rotation failed");
  }

  {
    Geometry::QuaternionRotateEigen qrot;
    a.clear();
    a.x()=1.;
    qrot.setAxis( geoCyl, Point(0,0,0), Point(0,1,0), pc::pi/2); // rotate around y-axis
    a = qrot(a); // rot. 90 deg.
    assert( eq(a.x(),0,1e-8) && "Vector rotation failed");
    a = qrot(a); // rot. 90 deg.
    assert( eq(a.x(),-1,1e-8) && "Vector rotation failed");
  }


  // check table of averages
  typedef Analysis::Table2D<float,Average<float> > Ttable;
  Ttable table(0.1, Ttable::XYDATA);
  table(2.1)+=1;
  table(2.1)+=3;
  assert( eq( table(2.1), 2 ) && "Bad average or 2D table");

  // check vector operations
  {
    Point a(2,3,5);
    assert( eq(a.len(), sqrt(4+9+25), 1e-8) && "Bad len() calculation");
    a*=2;
  }

}
