#include "faunus/geometry.h"
#include "faunus/inputfile.h"
#include "faunus/species.h"
#include "faunus/physconst.h"
#include "faunus/group.h"
#include "faunus/point.h"

namespace Faunus {
  namespace Geometry {

    double geometrybase::volume;

    double geometrybase::_dist(const point &p1, const point &p2) {
      return sqrt( _sqdist(p1, p2) );
    }

    void geometrybase::scale(point &a, const double &s) const {
    }

    bool geometrybase::collision(const particle &a, collisiontype type) {
      /*
         if (type==MATTER)
         for (int i=0; i<p.size(); i++)
         if (&p[i]!=&a) // avoid possible self overlap
         if ( a.overlap( p[i], _sqdist(a,p[i]) )==true )
         return true;
         */
      return false;
    }

    void geometrybase::pad(std::ostringstream& o, char w) {
      o << "#   " << setw(w) << std::left;
    }

    string geometrybase::info(char w) {
      std::ostringstream o;
      pad(o,w); o << "Boundary" << name << endl;
      return o.str();
    }

    void geometrybase::setvolume(double vol) {
      volume=vol;
    }

    double geometrybase::getvolume() const {
      return volume;
    }

    bool geometrybase::save(string file) {
      std::ofstream fout( file.c_str());
      if (fout) {
        fout.precision(10);
        fout <<  getvolume() << endl;
        return true;
      }
      return false;
    }

    /*!
     * \param file Filename
     * \param resize True if the current geometry should be resized to match file content (default: false)
     */
    bool geometrybase::load(string file, bool resize) {
      double vol;
      std::ifstream f( file.c_str() );
      if (f) {
        f >> vol;
        setvolume(vol);
        f.close();
        return true;
      }
      std::cerr << "# Geometry data NOT read from file " << file << endl;
      return false;
    }

    //
    //--- sphere geometry ---
    //

    void sphere::setvolume(double vol) {
      volume=vol;
      setradius( pow( 3*vol/(4*pc::pi), 1/3.) );
    }

    sphere::sphere(double radius) {
      setradius(radius);
    }

    sphere::sphere(inputfile &in)  {
      setradius(in.getflt("sphere_radius"));
    }

    void sphere::setradius(double radius) {
      name="Spherical";
      r = radius; 
      r2 = r*r; 
      diameter = 2*r; 
      volume = (4./3.)*pc::pi*r*r*r;
    }

    string sphere::info(char w) {
      std::ostringstream o;
      o << Geometry::geometrybase::info(w);
      pad(o,w); o << "Radius" << r << endl;
      return o.str();
    }

    void sphere::randompos(point &m) {
      double l=r2+1;
      while (l>r2) {
        m.x = slp.random_half()*diameter;
        m.y = slp.random_half()*diameter;
        m.z = slp.random_half()*diameter;
        l=m.x*m.x+m.y*m.y+m.z*m.z;
      }
    }

    point sphere::vdist(const point&a, const point&b) {
      return a-b;
    }


    void sphere::boundary(point &m) const {}

    bool sphere::collision(const particle &a, collisiontype type) {
      return (a.x*a.x+a.y*a.y+a.z*a.z > r2) ? true:false;
    }

    //
    //--- cuboid geometry ---
    //

    point cuboid::len;
    point cuboid::len_half;
    point cuboid::len_inv;

    cuboid::cuboid(inputfile &in) {
      name="Cuboid";
      double cubelen=in.getflt("cuboid_len",-1);
      if (cubelen<=0) {
        len.x=in.getflt("cuboid_xlen",0);
        len.y=in.getflt("cuboid_ylen",0);
        len.z=in.getflt("cuboid_zlen",0);
      } else {
        len.x=cubelen;
        len.y=cubelen;
        len.z=cubelen;
      }
      setlen(len);
      point min;
      min.x=in.getflt("cuboid_xmin",0);
      min.y=in.getflt("cuboid_ymin",0);
      min.z=in.getflt("cuboid_zmin",0);
      point max;
      max.x=in.getflt("cuboid_xmax",len.x);
      max.y=in.getflt("cuboid_ymax",len.y);
      max.z=in.getflt("cuboid_zmax",len.z);
      setslice(min,max);
    }

    bool cuboid::setlen(point l) {
      assert(l.x>0);              // debug information
      assert(l.y>0);              // 
      assert(l.z>0);              // 
      if (l.x<=0  ||              // check non-negative value
          l.y<=0  ||              // 
          l.z<=0  )               // 
        return false;             // 
      len = l;                    // cuboid sidelength
      len_half=l*0.5;             // half cuboid sidelength
      len_inv.x=1./len.x;         // inverse cuboid side length
      len_inv.y=1./len.y;         // 
      len_inv.z=1./len.z;         // 
      volume = len.x*len.y*len.z; // cuboid volume
      return true;
    }

    bool cuboid::setslice(point min, point max) {
      assert(min.x>=0);              // debug information
      assert(min.y>=0);              // 
      assert(min.z>=0);              // 
      if (min.x<0  ||                // check non-negative value
          min.y<0  ||                // 
          min.z<0  )                 // 
        return false;                // 
      assert(max.x<=len.x);          // debug information
      assert(max.y<=len.y);          // 
      assert(max.z<=len.z);          // 
      if (max.x>len.x  ||            // check non-negative value
          max.y>len.y  ||            // 
          max.z>len.z  )             // 
        return false;                // 
      slice_min = len_half-max;      // set minimum corner (other way around than in input!!)
      slice_max = len_half-min;      // set maximum corner
      return true;
    }

    string cuboid::info(char w) {
      std::ostringstream o;
      o << Geometry::geometrybase::info(w);
      pad(o,w); o << "Sidelengths" << len.x << " x " << len.y << " x " << len.z << endl;
      pad(o,w); o << "Slice position [x y z]"
        << len_half.x-slice_max.x << "-" << len_half.x-slice_min.x << " " 
        << len_half.y-slice_max.y << "-" << len_half.y-slice_min.y << " "
        << len_half.z-slice_max.z << "-" << len_half.z-slice_min.z << endl;
      return o.str();
    }

    point cuboid::randompos() {
      point m;
      randompos(m);
      return m;
    }

    void cuboid::randompos(point &m) {
      m.x = slp.random_half()*len.x;
      m.y = slp.random_half()*len.y;
      m.z = slp.random_half()*len.z;
    }

    bool cuboid::collision(const particle &a, collisiontype type) {
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

    bool cuboid::save(string file) {
      if ( Geometry::geometrybase::save(file) ) {
        std::ofstream fout( file.c_str(), std::ios_base::app);
        if (fout) {
          fout.precision(10);
          fout << len.x << " " << len.y << " " << len.z << endl;
          return true;
        }
      }
      return false;
    }

    bool cuboid::load(string file, bool resize) {
      point l;
      if ( Geometry::geometrybase::load(file, resize) ) {
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


    //
    //--- cuboid slit geometry ---
    //

    cuboidslit::cuboidslit(inputfile &in) : cuboid(in) {
      name="Cuboid XY-periodicity";
    }

    /*!
     * \param length Length of the cylinder
     * \param radius Radius of the cylinder
     */
    cylinder::cylinder(double length, double radius) {
      init(length, radius);
    }

    cylinder::cylinder(inputfile &in) {
      init( in.getflt("cylinder_len", 0), in.getflt("cylinder_radius", 0) );
    }

    void cylinder::init(double length, double radius) {
      name="Cylindrical";
      len=length;
      r=radius;
      r2=r*r;
      diameter=r*2;
      volume=2*r2*pc::pi*len;
      halflen=len/2;
    }

    void cylinder::randompos(point &m) {
      double l=r2+1;
      m.z = slp.random_half()*len;
      while (l>r2) {
        m.x = slp.random_half()*diameter;
        m.y = slp.random_half()*diameter;
        l=m.x*m.x+m.y*m.y;
      }
    }

    bool cylinder::collision(const particle &a, collisiontype type) {
      return 
        ( a.x*a.x+a.y*a.y>r2 || ( a.z<-halflen || a.z>halflen ) ) ? true:false;
    }

    string cylinder::info(char w) {
      std::ostringstream o;
      o << Geometry::geometrybase::info(w);
      pad(o,w); o << "Length (A)" << len << endl;
      pad(o,w); o << "Radius (A)" << r << endl;
      return o.str();
    }

#ifdef HYPERSPHERE
    const double hypersphere::pi=pc::pi;

    hypersphere::hypersphere(inputfile &in) : sphere(in) {
    }

    bool hypersphere::collision(const particle &p) {
      return false;
    }

    void hypersphere::randompos(point &m) {
      double rho=sqrt(slp.random_one());
      double omega=slp.random_one()*2.*pi;
      double fi=slp.random_one()*2.*pi;
      m.z1=sqrt(1.-rho*rho);
      m.z2=m.z1*cos(omega);
      m.z1=m.z1*sin(omega);
      m.z3=rho*sin(fi);
      m.z4=rho*cos(fi);
    }

    string hypersphere::info() {
      std::ostringstream o;
      o << Geometry::geometrybase::info() 
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
  }//namespace geometry

}//namespace
