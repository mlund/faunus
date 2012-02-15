#include "faunus/geometry.h"
#include "faunus/inputfile.h"
#include "faunus/species.h"
#include "faunus/physconst.h"
#include "faunus/group.h"
#include "faunus/point.h"
#include <faunus/faunus.h>

namespace Faunus {
  namespace Geometry {

    using namespace Faunus::textio;

    Geometrybase::~Geometrybase() {}

    double Geometrybase::dist(const Point &p1, const Point &p2) {
      return sqrt(sqdist(p1,p2));
    }

    void Geometrybase::scale(Point &a, const double &s) const {
      assert(!"Volume scaling function unimplemented for this geometry");
    }

    string Geometrybase::info(char w) {
      using namespace textio;
      std::ostringstream o;
      o << pad(SUB,w, "Boundary") << name << endl
        << pad(SUB,w, "Volume") << getVolume() << _angstrom << cubed
        << " = " << getVolume()/1e27 << " liters" << endl
        << _info(w);
      return o.str();
    }

    void Geometrybase::setVolume(double volume) {
      assert( volume>0 && "Zero volume not allowed!");
      _setVolume( volume );
      assert( std::abs( (volume-getVolume())/volume )<1e-9 && "setVolume() and/or getVolume() is broken!" );
    }

    double Geometrybase::getVolume() const {
      return _getVolume();
    }

    bool Geometrybase::save(string file) {
      std::ofstream fout( file.c_str());
      if (fout) {
        fout.precision(10);
        fout <<  getVolume() << endl;
        return true;
      }
      return false;
    }

    /*!
     * \param file Filename
     * \param resize True if the current geometry should be resized to match file content (default: false)
     */
    bool Geometrybase::load(string file, bool resize) {
      double vol;
      std::ifstream f( file.c_str() );
      if (f) {
        f >> vol;
        setVolume(vol);
        f.close();
        return true;
      }
      std::cerr << "!! Geometry data NOT read from file " << file << endl;
      return false;
    }

    Sphere::Sphere(double radius) {
      setradius(radius);
    }

    Sphere::Sphere(InputMap &in, string prefix)  {
      setradius( in.get<double>(prefix+"_radius", -1.0) );
    }

    void Sphere::setradius(double radius) {
      assert(radius>0 && "Radius must be larger than zero.");
      name="Spherical";
      r = radius; 
      r2 = r*r; 
      diameter = 2*r; 
    }

    double Sphere::_getVolume() const {
      return 4/3.*pc::pi*std::pow(r,3);
    }

    void Sphere::_setVolume(double vol) {
      setradius( pow( 3*vol/(4*pc::pi), 1/3.) );
    }

    void Sphere::scale(Point &a, const double &newvolume) const {
      assert( getVolume()>0 && newvolume>0 );
      double newradius = pow( 3*newvolume/(4*pc::pi), 1/3.);
      a = a * (newradius/r);
    }

    string Sphere::_info(char w) {
      std::ostringstream o;
      o << pad(SUB,w,"Radius") << r << textio::_angstrom << endl;
      return o.str();
    }

    void Sphere::randompos(Point &m) {
      double l=r2+1;
      while (l>r2) {
        m.x = slp.randHalf()*diameter;
        m.y = slp.randHalf()*diameter;
        m.z = slp.randHalf()*diameter;
        l=m.x*m.x+m.y*m.y+m.z*m.z;
      }
    }

    Point Sphere::vdist(const Point&a, const Point&b) {
      return a-b;
    }

    void Sphere::boundary(Point &m) const {}

    bool Sphere::collision(const particle &a, collisiontype type) const {
      return (a.x*a.x+a.y*a.y+a.z*a.z > r2) ? true:false;
    }

    //
    //--- Cuboid geometry ---
    //

    Cuboid::Cuboid(InputMap &in) {
      name="Cuboid";
      double cubelen=in.get<double>("cuboid_len",-1);
      if (cubelen<=0) {
        len.x=in.get<double>("cuboid_xlen",0);
        len.y=in.get<double>("cuboid_ylen",0);
        len.z=in.get<double>("cuboid_zlen",0);
      } else {
        len.x=cubelen;
        len.y=cubelen;
        len.z=cubelen;
      }
      setlen(len);
      Point min;
      min.x=in.get<double>("cuboid_xmin",0);
      min.y=in.get<double>("cuboid_ymin",0);
      min.z=in.get<double>("cuboid_zmin",0);
      Point max;
      max.x=in.get<double>("cuboid_xmax",len.x);
      max.y=in.get<double>("cuboid_ymax",len.y);
      max.z=in.get<double>("cuboid_zmax",len.z);
      setslice(min,max);
    }

    bool Cuboid::setlen(Point l) {
      assert(l.x>0 && l.y>0 && l.z>0);
      if (l.x<=0||l.y<=0||l.z<=0) 
        return false;
      len = l;                    // Cuboid sidelength
      len_half=l*0.5;             // half Cuboid sidelength
      len_inv.x=1/len.x;          // inverse Cuboid side length
      len_inv.y=1/len.y;
      len_inv.z=1/len.z;
      return true;
    }

    bool Cuboid::setslice(Point min, Point max) {
      assert(min.x>=0 && min.y>=0 && min.z>=0);
      assert(max.x<=len.x);          // debug information
      assert(max.y<=len.y);          // 
      assert(max.z<=len.z);          // 

      if (min.x<0  ||                // check non-negative value
          min.y<0  ||                // 
          min.z<0  )                 // 
        return false;                // 

      if (max.x>len.x  ||            // check non-negative value
          max.y>len.y  ||            // 
          max.z>len.z  )             // 
        return false;                // 
      slice_min = len_half-max;      // set minimum corner (other way around than in input!!)
      slice_max = len_half-min;      // set maximum corner
      return true;
    }

    void Cuboid::_setVolume(double newV) {
      scale(len,newV);
      setlen(len);
    }

    double Cuboid::_getVolume() const {
      return len.x*len.y*len.z;
    }

    string Cuboid::_info(char w) {
      std::ostringstream o;
      o << pad(SUB,w, "Sidelengths") << len.x << " x " << len.y << " x " << len.z << " ("+textio::angstrom+")" << endl
        << pad(SUB,w, "Slice position [x y z]")
        << len_half.x-slice_max.x << "-" << len_half.x-slice_min.x << " " 
        << len_half.y-slice_max.y << "-" << len_half.y-slice_min.y << " "
        << len_half.z-slice_max.z << "-" << len_half.z-slice_min.z << endl;
      return o.str();
    }

    Point Cuboid::randompos() {
      Point m;
      randompos(m);
      return m;
    }

    void Cuboid::randompos(Point &m) {
      m.x = slp.randHalf()*len.x;
      m.y = slp.randHalf()*len.y;
      m.z = slp.randHalf()*len.z;
    }

    bool Cuboid::collision(const particle &a, collisiontype type) const {
      switch (type) {
        case (BOUNDARY): // collision with geometry boundaries
          if (std::abs(a.x) > len_half.x ||
              std::abs(a.y) > len_half.y ||
              std::abs(a.z) > len_half.z  ) return true;
          break;
        case (ZONE):     // collision with forbidden zone (slice)
          if (std::abs(a.x) > len_half.x  ||
              std::abs(a.y) > len_half.y  ||
              std::abs(a.z) > len_half.z  ||
              a.x  < slice_min.x ||
              a.y  < slice_min.y ||
              a.z  < slice_min.z ||
              a.x  > slice_max.x ||
              a.y  > slice_max.y ||
              a.z  > slice_max.z  )
          return true;
          break;
      }
      return false;
    }

    bool Cuboid::save(string file) {
      if ( Geometrybase::save(file) ) {
        std::ofstream fout( file.c_str(), std::ios_base::app);
        if (fout) {
          fout.precision(10);
          fout << len.x << " " << len.y << " " << len.z << endl;
          return true;
        }
      }
      return false;
    }

    bool Cuboid::load(string file, bool resize) {
      Point l;
      if ( Geometrybase::load(file, resize) ) {
        std::ifstream f( file.c_str() );
        if (f) {
          f >> l.x >> l.y >> l.z;
          setlen(l);
          f.close();
          return true;
        }
      }
      std::cerr << "# Container data NOT read from file " << file << endl;
      return false;
    }

    void Cuboid::scale(Point &a, const double &newvolume) const {
      assert( getVolume()>0 && newvolume>0 );
      a = a * std::pow( newvolume/getVolume(), 1/3.);
    }

    Cuboidslit::Cuboidslit(InputMap &in) : Cuboid(in) {
      name="Cuboid XY-periodicity";
    }

    /*!
     * \param length Length of the Cylinder (angstrom)
     * \param radius Radius of the Cylinder (angstrom)
     */
    Cylinder::Cylinder(double length, double radius) {
      init(length, radius);
    }

    /*!
     * The inputmap is searched for the keywords:
     * \li \c cylinder_len (along z) [A]
     * \li \c cylinder_radius [A]
     */
    Cylinder::Cylinder(InputMap &in) {
      init( in.get<double>("cylinder_len", 0), in.get<double>("cylinder_radius", 0) );
    }

    void Cylinder::init(double length, double radius) {
      name="Cylindrical (hard ends)";
      assert(length>0 && radius>0 && "Cylinder length and radius must be bigger than zero.");
      len=length;
      setVolume( pc::pi*radius*radius*len );
    }

    void Cylinder::_setVolume(double newV) {
      r2=newV/(pc::pi*len);
      r=sqrt(r2);
      diameter=2*r;
      halflen=len/2;
    }

    double Cylinder::_getVolume() const {
      return r2*pc::pi*len;
    }

    void Cylinder::boundary(Point &p ) const {}

    void Cylinder::randompos(Point &m) {
      double l=r2+1;
      m.z = slp.randHalf()*len;
      while (l>r2) {
        m.x = slp.randHalf()*diameter;
        m.y = slp.randHalf()*diameter;
        l=m.x*m.x+m.y*m.y;
      }
    }

    bool Cylinder::collision(const particle &a, collisiontype type) const {
      assert( (halflen-len/2)<1e-9 && "Cylinder length initialization problems" );
      return 
        ( a.x*a.x+a.y*a.y>r2 || ( a.z<-halflen || a.z>halflen ) ) ? true:false;
    }

    string Cylinder::_info(char w) {
      std::ostringstream o;
      o << pad(SUB,w, "Length") << len << textio::_angstrom << endl
        << pad(SUB,w, "Radius") << r << textio::_angstrom << endl;
      return o.str();
    }

    PeriodicCylinder::PeriodicCylinder(double length, double radius) : Cylinder(length,radius) {
      name="Cylindrical (periodic ends)";
    }

    PeriodicCylinder::PeriodicCylinder(InputMap &in) : Cylinder(in) {
      name="Cylindrical (periodic ends)";
    }

    void PeriodicCylinder::boundary(Point &a) const {
      if (std::abs(a.z)>halflen)
        a.z-=len*anint(a.z/len);
    }

#ifdef HYPERSPHERE
    const double hyperSphere::pi=pc::pi;

    hyperSphere::hyperSphere(InputMap &in) : Sphere(in) {
    }

    bool hyperSphere::collision(const particle &p) const {
      return false;
    }

    void hyperSphere::randompos(point &m) {
      double rho=sqrt(slp.random_one());
      double omega=slp.random_one()*2.*pi;
      double fi=slp.random_one()*2.*pi;
      m.z1=sqrt(1.-rho*rho);
      m.z2=m.z1*cos(omega);
      m.z1=m.z1*sin(omega);
      m.z3=rho*sin(fi);
      m.z4=rho*cos(fi);
    }

    string hyperSphere::_info() {
      std::ostringstream o;
      o << "#   Shape                = Hyperspherical" << endl
        << "#   Radius               = " << r << endl;
      return o.str();
    }

    /*
       void hyperSphere::move(int i, double dangle) {
       double nfi=2.*acos(-1.)*slp.random_one();
       double nrho=sqrt((slp.random_one()-1.)*sin(dangle)*sin(dangle)+1.);
       double nomega=(2.*slp.random_one()-1.)*acos(cos(dangle)/nrho);
       move(trial[i],nrho,nomega,nfi);
       }
       */
#endif

    Point massCenter(const Geometrybase &geo, const p_vec &p) {
      double sum=0;
      Point cm,t,o = p.at( p.size()/2 );  // set origo to middle particle
      for (auto &pi : p) {
        t = pi-o;              // translate to origo
        geo.boundary(t);        // periodic boundary (if any)
        cm += t * pi.mw;
        sum += pi.mw;
      }
      assert(sum>0 && "Zero mass no good in CM calculation! Did you forget to assign atom weights?");
      cm=cm*(1/sum) + o;
      geo.boundary(cm);
      return cm;
    }

    void translate(const Geometrybase &geo, p_vec &p, Point d) {
      for (auto &pi : p) {
        pi += d;
        geo.boundary(pi);
      }
    }

    FindSpace::FindSpace() {
      dir.x=1;
      dir.y=1;
      dir.z=1;
      allowContainerOverlap=false;
      allowMatterOverlap=false;   
    }

    bool FindSpace::matterOverlap(const Geometrybase &geo, const p_vec &p1, const p_vec &p2) {
      if (allowMatterOverlap==false)
        for (auto &i : p1)
          for (auto &j : p2) {
            double max=i.radius+j.radius;
            if ( geo.sqdist(i,j)<max*max )
              return true;
          }
      return false;
    }

    bool FindSpace::containerOverlap(const Geometrybase &geo, const p_vec &p) {
      if (allowContainerOverlap==false)
        for (auto &i : p)
          if (geo.collision(i)) return true;
      return false;
    }

    /*!
     * \param geo Geometry to use
     * \param dst Destination particle vector (will not be touched!)
     * \param p Particle vector to find space for. Coordinates will be changed.
     * \param maxtrials Number of times to try before timeout.
     */
    bool FindSpace::find(Geometrybase &geo, const p_vec &dst, p_vec &p, unsigned int maxtrials) {
      using namespace textio;
      cout << "Trying to insert " << p.size() << " particle(s)";
      Point cm,v;
      do {
        cout << ".";
        maxtrials--;
        cm = massCenter(geo, p);
        geo.randompos(v);
        v.x*=dir.x;
        v.y*=dir.y;
        v.z*=dir.z;
        translate(geo, p, -cm+v);
      } while (maxtrials>0 && (containerOverlap(geo,p)==true || matterOverlap(geo,p,dst)==true));
      if (maxtrials>0) {
        cout << " OK!\n";
        return true;
      }
      cout << " timeout!\n";
      assert(!"Timeout - found no space for particle(s).");
      return false;
    }

    /*!
     * \param geo Simulation geometry
     * \param beg Starting point of line to rotate around - typically a molecular mass center
     * \param end Ending point of line to rotate around
     * \param angle How many degrees to rotate
     */
    void VectorRotate::setAxis(Geometrybase &geo, const Point &beg, const Point &end, double angle) {
      assert(&geo!=nullptr);
      origin=beg;
      u=end-beg;
      geo.boundary(u);
      u=u*(1/geo.dist(beg,end));
      cosang=cos(angle);
      sinang=sin(angle);
      e1mcox=(1.-cosang)*u.x;
      e1mcoy=(1.-cosang)*u.y;
      e1mcoz=(1.-cosang)*u.z;
    }

    /*!
     * Rotate point around axis specified above.
     * \param geo Geometry
     * \param p Point to rotate
     */
    Point VectorRotate::rotate(const Geometrybase &geo, Point p) const {
      assert(&geo!=nullptr);
      Point b=p-origin;
      geo.boundary(b);           // Apply boundary conditions
      double eb=u.x*b.x + u.y*b.y + u.z*b.z;
      p.x=e1mcox*eb+cosang*b.x+sinang*(u.y*b.z - u.z*b.y) + origin.x;
      p.y=e1mcoy*eb+cosang*b.y+sinang*(u.z*b.x - u.x*b.z) + origin.y;
      p.z=e1mcoz*eb+cosang*b.z+sinang*(u.x*b.y - u.y*b.x) + origin.z;
      geo.boundary(p);
      return p;
    }
  }//namespace geometry

}//namespace
