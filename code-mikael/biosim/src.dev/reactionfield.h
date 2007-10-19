#include "point.h"
#include "interact.h"
#include "legendre.h"

class rfield {
 private:
  int steps;          // # of steps in summation (precision)
  double costheta;// cos(theta)
  legendre l;     // legendre polynomium class constructor

 public:
  double a;       // cavity radius
  double eo, ei;  // dielectric constants (solvent, inside sphere)
  rfield(double radius, double eps_o, double eps_i, int n=50) {
    a=radius;
    eo=eps_o;
    ei=eps_i;
    steps=n;  //number of series. Higher->better/slower.
    l.resize(steps);
  };

  void info() {
    cout << "  Dielectric constant (in,out) = " << ei << " " << eo << endl
         << "  Sphere radius                = " << a << endl
         << "  Summation length             = " << steps << endl;
  };

  //calc. potential in point "p" from charge in "p0".
  //ref: Woodward+Svensson, J.Phys.Chem, 1991, 95, p7471.
  //Note: multiply w. lB*eps_s to get kT.
  inline double phi(particle &p0, point &p, bool self=false) {
    if (p0.charge==0)
      return 0;
    double sum = 0;
    double r   = p.len();
    double r0  = p0.len();
    double dn;
    costheta=p.dot(p0) / (r * r0);
    l.eval(costheta); // legendre poly.

    //charge outside, potential outside
    //checked: OK!
    if (r0>a && r>a) {
      for (int n=1; n<steps; n+=1) {
	dn   = double(n);
	sum += pow(a*a/(r*r0),n+1) * l.p[n]
	  / (  (ei+eo*(1+1/dn))   );
      };
      sum *= (eo-ei) / (eo * a );
      if (self==false)
        sum += 1/(eo*p.dist(p0));
    };

    //charge inside, potential inside
    //checked: OK!
    if (r0<a && r<a) {
      for (int n=0; n<steps; n+=1) {
	dn = double(n);
	sum += pow(r*r0/(a*a),n) * l.p[n]
	  / ( eo-ei*(dn/(dn+1)) )  ;
      };
      sum *= (ei-eo) / (ei*a);
      if (self==false)
        sum += 1/(ei*p.dist(p0));
    };

    //charge inside, potential outside
    //checked: OK!
    if (r0<a && r>a) {
      for (int n=0; n<steps; n+=1) {
	dn=double(n);
	sum += (2*dn+1) * pow(r0/r,n) * l.p[n]
	  / (  ( dn*ei + eo*(dn+1) )  );
      };
      sum = sum / r;
    };

    //charge outside, potential inside (as above...)
    //checked: OK!
    if (r0>a && r<a) {
      for (int n=1; n<steps; n+=1) {
	dn = double(n);
	sum += pow(r/r0,n) * l.p[n]
	  / ( (ei+eo*(1+1/dn))   ) ;
      };
      sum *= (eo-ei) / (eo * r0) ;
      sum += 1/(eo*p.dist(p0));
    };
    
    return sum * p0.charge;
  };

};

class interact_rf : public rfield, public interact {
public:
  interact_rf(double bjerrum, double radius, double eps_s, double eps_i)
    : interact(bjerrum), rfield(radius,eps_s,eps_i) {};

  // born energy (vacuum->dielectric)
  double born(particle &p) {
    double eps=eo;
    if (p.len()<a)
      eps=ei;
    return -(lB*eo) * p.charge * p.charge / (2*p.radius) * (1-1./eps);
  };

  
  // self energy. "image charge"
  inline double selfenergy(particle &p) {
     return 0.5 * p.charge * phi(p,p,true);
  };

  inline double energy_rf(particle &p1, particle &p2) {
    return p2.charge * phi(p1,p2);
  };

  // energy of j'th particle. Units of kT.
  double energy_rf(vector<particle> &p, int j) {
    if (p[j].charge==0) return 0;
    int ps=p.size();
    double sum=0;
    for (int i=0; i<j; i++)
      sum += phi(p[i],p[j]);
    for (int i=j+1; i<ps; i++)
      sum += phi(p[i],p[j]);

    return lB * eo * ( sum*p[j].charge + selfenergy(p[j]) );
    //return lB * eo * sum * p[j].charge;
  };

  // energy of ghost particle (kT) incl. self energy.
  double energy_rf(vector<particle> &p, particle &a) {
    if (a.charge==0) return 0;
    double sum=0;
    for (int i=0; i<p.size(); i++)
      sum += phi( p[i], a);
    return lB * eo * ( sum*a.charge + selfenergy(a) );
  };

  // Total electrostatic potential in a point
  double potential_rf(vector<particle> &p, point &a) {
    double x=0;
    int ps=p.size();  
    for (int i=0; i<ps; i++)
      x+=phi( p[i], a );
    return lB*eo*x;
  };

};

