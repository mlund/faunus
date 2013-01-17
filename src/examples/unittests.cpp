#undef NDEBUG
#include <faunus/faunus.h>

using namespace Faunus;

bool eq(double a, double b, double tol=1e-6) { return (std::abs(a-b)<tol) ? true : false; }

/*
 * Check various copying operations
 * between particle types
 */
template<typename Tparticle>
void checkParticle() {
  Tparticle p;
  p.clear();
  p.mw=1.0; // set some property - should should not be overridden!
  Point::Tvec vec(1,0,0);
  p=vec;
  assert( eq(p.norm(), 1) );
  p+=vec;
  assert( eq(p.norm(), 2) );
  p*=0.5;
  assert( eq(p.norm(), 1) );

  assert( eq(p.mw, 1) && "mw should not be overridden!");

  p=Point::Tvec(2,3,5);
  assert( eq(p.len(), sqrt(4+9+25), 1e-8) && "Bad len() calculation");

  assert( eq(p.mw, 1) && "mw should not be overridden!");

  PointParticle a;
  a.clear();
  a.x()=10.;
  a.charge=1.0;
  p=a; 
  assert( eq(p.x(), a.x()) );
  assert( eq(p.charge, a.charge) );
  assert( eq(p.mw, 0) );
}

int main() {

  // check particle operations
  {
    checkParticle<PointParticle>();
    checkParticle<DipoleParticle>();
    checkParticle<CigarParticle>();
#ifdef FAU_HYPERSPHERE
    checkParticle<HyperParticle>();
#endif
  }

  // check infinity
  assert( pc::infty == -std::log(0) && "Problem with infinity" );

  // check approximate exp()
  assert( eq( exp_cawley(1),std::exp(1),1e-1 ) && "Problem with approximate exp() function" );

  // check approximate 1/sqrt()
  assert( eq( invsqrtQuake(20.), 1/std::sqrt(20.), 1e-1) && "Problem w. inverse Quake sqrt.");

  // groups
  {
    Group g(2,5);           // first, last particle                                                                  
    assert(g.front()==2);
    assert(g.back()==5);
    assert(g.size()==4);

    g.resize(0);
    assert(g.empty());

    g.resize(1000);
    assert(g.size()==1000);

    g.setfront(1);
    g.setback(10);
    assert(g.front()==1);
    assert(g.back()==10);
    assert(g.size()==10);

    int cnt=0;
    for (auto i : g) {
      cnt++;
      assert(g.find(i));
    }
    assert(cnt==g.size());

    assert(!g.find(0));
    assert(!g.find(11));

    // check random
    int min=1e6, max=-1e6;
    for (int n=0; n<1e6; n++) {
      int i=g.random();
      if (i<min) min=i;
      if (i>max) max=i;
      assert( g.find(i) );
    }
    assert(min=g.front());
    assert(max=g.back());
  }

  // check geometries
  {
    Geometry::Sphere geoSph(1000);
    Geometry::Cylinder geoCyl(1000,1000);
    Point a( 0,0,sqrt(16)), b(0,sqrt(64),0);
    double x = geoSph.sqdist(a,b);
    double y = geoCyl.sqdist(a,b);
    assert( x==16+64 && y==16+64 );
    assert( eq(x,y) && "No good distance calculation");
  }

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
    Geometry::Cylinder geo(1000,1000);
    Point a(0,0,0);
    a.x()=1.;
    vrot.setAxis( geo, Point(0,0,0), Point(0,1,0), pc::pi/2); // rotate around y-axis
    a = vrot(a); // rot. 90 deg.
    assert( eq(a.x(),0,1e-8) && "Vector rotation failed");
    a = vrot(a); // rot. 90 deg.
    assert( eq(a.x(),-1,1e-8) && "Vector rotation failed");
  }

  {
    Geometry::QuaternionRotate qrot;
    Geometry::Cylinder geo(1000,1000);
    Point a(0,0,0);
    a.x()=1.;
    qrot.setAxis( geo, Point(0,0,0), Point(0,1,0), pc::pi/2); // rotate around y-axis
    a = qrot(a); // rot. 90 deg.
    assert( eq(a.x(),0,1e-8) && "Vector rotation failed");
    a = qrot(a); // rot. 90 deg.
    assert( eq(a.x(),-1,1e-8) && "Vector rotation failed");
  }

  {
    Geometry::QuaternionRotateEigen qrot;
    Geometry::Cylinder geo(1000,1000);
    Point a(0,0,0);
    a.x()=1.;
    qrot.setAxis( geo, Point(0,0,0), Point(0,1,0), pc::pi/2); // rotate around y-axis
    a = qrot(a); // rot. 90 deg.
    assert( eq(a.x(),0,1e-8) && "Vector rotation failed");
    a = qrot(a); // rot. 90 deg.
    assert( eq(a.x(),-1,1e-8) && "Vector rotation failed");
  }

  // check table of averages
  {
    typedef Analysis::Table2D<float,Average<float> > Ttable;
    Ttable table(0.1, Ttable::XYDATA);
    table(2.1)+=1;
    table(2.1)+=3;
    assert( eq( table(2.1), 2 ) && "Bad average or 2D table");
  }

}
