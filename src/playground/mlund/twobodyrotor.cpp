#include <faunus/faunus.h>

using namespace Faunus;
using namespace Faunus::Potential;
typedef Space<Geometry::Sphere> Tspace;
typedef CombinedPairPotential<DebyeHuckel,LennardJonesLB> Tpairpot;

// Class for spherical coordinates
// more info:
// http://blog.demofox.org/2013/10/12/converting-to-and-from-polar-spherical-coordinates-made-easy/
struct Spherical {

  double r,theta,phi;

  Spherical() : r(0), theta(0), phi(0) {}

  Spherical(double radius, double polar=0, double azimuth=0) {
    r=radius;
    theta=polar;
    phi=azimuth;
  }

  Spherical& operator=(const Point &c) {
    r=c.norm();
    theta=atan2(c.y(),c.x());
    phi=acos(c.z()/r);
    return *this;
  }

  Point cartesian() const {
    return r*Point( cos(theta)*cos(phi), sin(theta)*cos(phi), sin(phi) );
  }
};

int main(int argc, char** argv) {
  InputMap mcp("twobodyrotor.input");
  Tspace spc(mcp);
  spc.geo.setRadius(1e5);
  auto pot = Energy::Nonbonded<Tspace,Tpairpot>(mcp);

  double rmin = mcp("rmin",0.);
  double rmax = mcp("rmax",50.);
  double dr   = mcp("dr",0.5);
  double dangle = mcp("dangle",0.01);

  double Q=0; // partition function
  double U=0; // energy sum

  long unsigned int cnt=0;

  // load two molecules
  vector<Group> pol(2);
  Tspace::ParticleVector v;
  FormatAAM::load("mol.aam", v);
  Geometry::cm2origo(spc.geo, v);
  for (auto &i : pol) {
    i = spc.insert(v);
    spc.enroll(i);
    i.setMassCenter(spc);
  }

  // loop over distances
  for (double r=rmin; r<rmax; r+=dr) {
    pol[1].translate(spc, -pol[1].cm + Point(0,0,r));
    pol[1].accept(spc);

    // scan angles...
    for (double x1=0; x1<pc::pi; x1+=dangle) {

      Point a = Spherical(1,x1).cartesian();

      for (double x2=0; x2<2*pc::pi; x2+=dangle) {
        pol[0].rotate(spc, a, dangle);
        pol[0].accept(spc);

        for (double y1=0; y1<pc::pi; y1+=dangle) {

          Point b = Spherical(1,y1).cartesian();

          for (double y2=0; y2<2*pc::pi; y2+=dangle) {
            pol[1].rotate(spc, b, dangle);
            pol[1].accept(spc);

            double betau = pot.g2g(spc.p, pol[0], pol[1]);
            Q += exp(-betau);
            U += betau * exp(-betau);
            cnt++;
          }
        }
      }
    }
    cout << r << " " << Q/cnt << " " << U/Q << "\n";
  } // end of distance scan

}
