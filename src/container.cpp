#include "faunus/container.h"

namespace Faunus {
  //----------- CONTAINER -----------------
  string container::info() {
    double z=charge();
    std::ostringstream o;
    o << endl
      << "# SIMULATION CONTAINER:" << endl
      << "#   Number of particles  = " << p.size() << endl
      << "#   Volume (AA^3)        = " << volume << endl
      << "#   Electroneutrality    = " 
      << ((abs(z)>1e-7) ? "NO!" : "Yes") << " "  << z << endl;
    return o.str();
  }
  string container::povray() {
    return string(
        "#declare cell=texture {\n"
        "pigment {color rgbf <1,1,1,.9>}\n"
        " finish {phong .1 ambient .2}\n"
        "}\n" );
  }

  //----------- CELL ----------------------
  cell::cell(double radius) { setradius(radius); }
  cell::cell(inputfile &in)  {
    atom.load(in);
    setradius(in.getflt("cellradius"));
  }
  void cell::setradius(double radius) {
    r = radius; 
    r2 = r*r; 
    diameter = 2*r; 
    volume = (4./3.)*acos(-1.)*r*r*r;
  }
  string cell::info() {
    std::ostringstream o;
    o << container::info() 
      << "#   Shape                = Spherical" << endl
      << "#   Radius               = " << r << endl;
    return o.str();
  }
  void cell::randompos(point &p) {
    double l=r2+1;
    while (l>r2) {
      p.x = slp.random_half()*diameter;
      p.y = slp.random_half()*diameter;
      p.z = slp.random_half()*diameter;
      l=p.x*p.x+p.y*p.y+p.z*p.z;
    }
  }
  string cell::povray() {
    std::ostringstream o;
    o << "sphere {<0,0,0>," << r << " texture {cell}}\n"
      << "cylinder {<0,0," << -r << ">,<0,0," << r << ">,0.5 texture {cell}}\n";
    return o.str();
  }
  //----------- BOX --------------------------
  void box::setvolume(double v) {
    len=pow(v, 1./3);
    setlen(len);
  }
  bool box::setlen(double l) {
    if (l<=0)
      return false;
    len = l;    // cubic box sidelength
    len_half=len/2.;
    len_inv=1./len;
    volume = len*len*len;
    return true;
  }
  box::box(double l) { setlen(l); }

  box::box(inputfile &in) {
    atom.load(in);
    if (!setlen( in.getflt("boxlen",-1) ))
      setvolume( in.getflt("volume") );
  }

  string box::info() {
    std::ostringstream o;
    o << container::info() 
      << "#   Shape                = Cube" << endl
      << "#   Side length          = " << len << endl;
    return o.str();
  }
  point box::randompos() {
    point p;
    randompos(p);
    return p;
  }
  void box::randompos(point &p) {
    p.x = slp.random_half()*len;
    p.y = slp.random_half()*len;
    p.z = slp.random_half()*len;
  }
  /*void box::randompos(vector<point> &p) {
    short unsigned size=p.size;
    }*/
  string box::povray() {
    std::ostringstream o;
    o << "box {<" <<-len_half <<"," <<-len_half <<"," <<-len_half <<"> , <"
      << len_half <<"," <<len_half <<"," <<len_half <<"> texture {cell}}\n";
    return o.str();
  }
  //----------- XYPLANE ---------------------
  xyplane::xyplane(inputfile &in) : box(in) {
  }
  void xyplane::randompos(point &p) {
    p.x = slp.random_half()*len;
    p.y = slp.random_half()*len;
    p.z = 0;
  }

  //----------- SLIT --------------------------
  slit::slit(inputfile &in) : box(in) {
    if (in.getflt("zboxlen", 0)>0)
      zlen=in.getflt("zboxlen");
    else
      zlen=in.getflt("boxlen");
    zlen_half=zlen*0.5;
    xyarea=pow(len, 2.);
    volume=xyarea*zlen;
  };
  string slit::info() {
    std::ostringstream o;
    o << container::info() 
      << "#   Shape                = Cube - xy periodicity, only" << endl
      << "#   XY-Side length          = " << len << endl
      << "#   Z -Side length          = " << zlen << endl;
    return o.str();
  }

  //-------------- CLUTCH -------------------------
  //! \param radius Radius of the cell
  //! \param min Beginning of the particle-excluded region
  //! \param max Ending of the particle-excluded region (Note max>min)
  clutch::clutch(double radius, double min, double max) {
    r=radius;
    r2=r*r;
    diameter=2*r;
    volume=(4./3.)*acos(-1.)*r2*r;
    zmin=min;
    zmax=max;
  }
  void clutch::randompos(point &p) {
    double l=r2+1;
    while (l>r2) {
      p.x = slp.random_half()*diameter;
      p.y = slp.random_half()*diameter;
      p.z = slp.random_half()*diameter;
      if (p.z>zmax || p.z<zmin)
        l=p.x*p.x+p.y*p.y+p.z*p.z; //squared distance from origo
    };
  }
  //------------CYLINDER---------------------------
  //! \param length    Length of the cylinder
  //! \param radius    Radius of the cylinder
  cylinder::cylinder(double length, double radius) {
    len=length;
    r=radius;
    r2=r*r;
    diameter=r*2;
    volume=2*r2*acos(-1.)*len;
  }
  void cylinder::randompos(point &p) {
    double l=r2+1;
    p.z = slp.random_one()*len;
    while (l>r2) {
      p.x = slp.random_half()*diameter;
      p.y = slp.random_half()*diameter;
      l=p.x*p.x+p.y*p.y;
    };
  }
  string cylinder::info() {
    std::ostringstream o;
    o << container::info()
      << "#   Shape                = Cylindrical" << endl
      << "#   Lenght               = " << len <<endl
      << "#   Radius               = " << r << endl;
    return o.str();
  }
  string cylinder::povray() {
    std::ostringstream o;
    o << "cylinder {<0,0,0>,<0,0" << len << ">," << r <<" texture {cell}}\n"
      << "cylinder {<0,0,0>,<0,0" << len << ">,0.5 texture {cell}}\n";
    return o.str();
  }     

#ifdef HYPERSPHERE
  const double hypersphere::pi=3.141592654;
  hypersphere::hypersphere(inputfile &in) : cell(in) {
  }
  bool hypersphere::collision(const particle &p) {
    return false;
  }
  void hypersphere::randompos(point &p) {
    double rho=sqrt(slp.random_one());
    double omega=slp.random_one()*2.*pi;
    double fi=slp.random_one()*2.*pi;
    p.z1=sqrt(1.-rho*rho);
    p.z2=p.z1*cos(omega);
    p.z1=p.z1*sin(omega);
    p.z3=rho*sin(fi);
    p.z4=rho*cos(fi);
  }
  string hypersphere::info() {
    std::ostringstream o;
    o << container::info() 
      << "#   Shape                = Hyperspherical" << endl
      << "#   Radius               = " << r << endl;
    return o.str();
  }
  /*
  void hypersphere::move(int i, double dangle) {
    double nfi=2.*acos(-1.)*slp.random_one();
    double nrho=sqrt((slp.random_one()-1.)*sin(dangle)*sin(dangle)+1.);
    double nomega=(2.*slp.random_one()-1.)*acos(cos(dangle)/nrho);
    move(trial[i],nrho,nomega,nfi);
  }
  */
#endif

}//namespace
