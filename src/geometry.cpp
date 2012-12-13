#include <faunus/geometry.h>
#include <faunus/inputfile.h>
#include <faunus/species.h>
#include <faunus/physconst.h>
#include <faunus/textio.h>
#include <faunus/group.h>

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
        << " = " << getVolume()/1e3 << " nm" << cubed
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
      std::ifstream f( file.c_str() );
      if (f) {
        double vol;
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
      setradius( in.get<double>(prefix+"_radius", -1.0, "Spherical container radius (A)") );
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
        m.x() = slp.randHalf()*diameter;
        m.y() = slp.randHalf()*diameter;
        m.z() = slp.randHalf()*diameter;
        l=m.x()*m.x()+m.y()*m.y()+m.z()*m.z();
      }
    }

    void Sphere::boundary(Point &m) const {}

    bool Sphere::collision(const particle &a, collisiontype type) const {
      return (a.x()*a.x()+a.y()*a.y()+a.z()*a.z() > r2) ? true:false;
    }

    //
    //--- Cuboid geometry ---
    //

    Cuboid::Cuboid(InputMap &in) {
      name="Cuboid";
      double cubelen=in.get<double>("cuboid_len",-1, name+" sidelength (AA)");
      if (cubelen<=0) {
        len.x()=in.get<double>("cuboid_xlen",0);
        len.y()=in.get<double>("cuboid_ylen",0);
        len.z()=in.get<double>("cuboid_zlen",0);
      } else {
        len.x()=cubelen;
        len.y()=cubelen;
        len.z()=cubelen;
      }
      setlen(len);
    }

    bool Cuboid::setlen(Point l) {
      assert(l.x()>0 && l.y()>0 && l.z()>0);
      if (l.x()<=0||l.y()<=0||l.z()<=0) 
        return false;
      len = l;                    // Cuboid sidelength
      len_half=l*0.5;             // half Cuboid sidelength
      len_inv.x()=1/len.x();          // inverse Cuboid side length
      len_inv.y()=1/len.y();
      len_inv.z()=1/len.z();
      return true;
    }

    void Cuboid::_setVolume(double newV) {
      scale(len,newV);
      setlen(len);
    }

    double Cuboid::_getVolume() const {
      return len.x()*len.y()*len.z();
    }

    string Cuboid::_info(char w) {
      std::ostringstream o;
      o << pad(SUB,w, "Sidelengths")
        << len.x() << " x " << len.y() << " x " << len.z() << " ("+textio::angstrom+")" << endl;
      return o.str();
    }

    Point Cuboid::randompos() {
      Point m;
      randompos(m);
      return m;
    }

    void Cuboid::randompos(Point &m) {
      m.x() = slp.randHalf()*len.x();
      m.y() = slp.randHalf()*len.y();
      m.z() = slp.randHalf()*len.z();
    }

    bool Cuboid::collision(const particle &a, collisiontype type) const {
      if (std::abs(a.x())>len_half.x()) return true;
      if (std::abs(a.y())>len_half.y()) return true;
      if (std::abs(a.z())>len_half.z()) return true;
      return false;
    }

    bool Cuboid::save(string file) {
      if ( Geometrybase::save(file) ) {
        std::ofstream fout( file.c_str(), std::ios_base::app);
        if (fout) {
          fout.precision(10);
          fout << len.x() << " " << len.y() << " " << len.z() << endl;
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
          f >> l.x() >> l.y() >> l.z();
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

    void Cylinder::boundary(Point &p) const {}

    void Cylinder::randompos(Point &m) {
      double l=r2+1;
      m.z() = slp.randHalf()*len;
      while (l>r2) {
        m.x() = slp.randHalf()*diameter;
        m.y() = slp.randHalf()*diameter;
        l=m.x()*m.x()+m.y()*m.y();
      }
    }

    bool Cylinder::collision(const particle &a, collisiontype type) const {
      assert( (halflen-len/2)<1e-9 && "Cylinder length initialization problems" );
      if ( a.z()<-halflen ) return true;
      if ( a.z()>halflen ) return true;
      if ( a.x()*a.x()+a.y()*a.y()>r2 ) return true;
      return false;
    }

    string Cylinder::_info(char w) {
      std::ostringstream o;
      o << pad(SUB,w, "Length") << halflen*2 << textio::_angstrom << endl
        << pad(SUB,w, "Radius") << sqrt(r2) << textio::_angstrom << endl;
      return o.str();
    }

    PeriodicCylinder::PeriodicCylinder(double length, double radius) : Cylinder(length,radius) {
      name="Cylindrical (periodic ends)";
    }

    PeriodicCylinder::PeriodicCylinder(InputMap &in) : Cylinder(in) {
      name="Cylindrical (periodic ends)";
    }

    void PeriodicCylinder::boundary(Point &a) const {
      if (std::abs(a.z())>halflen)
        a.z()-=len*anint(a.z()/len);
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
      if (p.empty())
        return Point();
      Group g(0, p.size()-1);
      return massCenter(geo,p,g);
    }

    Point massCenter(const Geometrybase &geo, const p_vec &p, const Group &g) {
      if (g.empty())
        return Point();
      assert(!p.empty());
      assert(g.back() < (int)p.size());
      assert(&geo!=NULL);
      double sum=0;
      Point cm(0,0,0);
      Point t;
      Point o = p.at( g.front()+(g.back()-g.front())/2 );  // set origo to middle particle
      for (auto i : g) {
        t = p[i]-o;              // translate to origo
        geo.boundary(t);        // periodic boundary (if any)
        cm += t * p[i].mw;
        sum += p[i].mw;
      }
      if (sum<1e-6) sum=1;
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

    void cm2origo(const Geometrybase &geo, p_vec &p) {
      translate(geo, p, -massCenter(geo, p) );
    }

    FindSpace::FindSpace() {
      dir.x()=1;
      dir.y()=1;
      dir.z()=1;
      allowContainerOverlap=false;
      allowMatterOverlap=false;   
    }

    FindSpace::~FindSpace() {}

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
        v.x()*=dir.x();
        v.y()*=dir.y();
        v.z()*=dir.z();
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
    
    VectorRotate::~VectorRotate() {}

    /*!
     * \param geo Simulation geometry
     * \param beg Starting point of line to rotate around - typically a molecular mass center
     * \param end Ending point of line to rotate around
     * \param angle Rotation angle [rad]
     */
    void VectorRotate::setAxis(Geometrybase &geo, const Point &beg, const Point &end, double angle) {
      assert(&geo!=nullptr);
      geoPtr=&geo;
      origin=beg;
      u=end-beg;
      geo.boundary(u);
      u=u*(1/geo.dist(beg,end));
      cosang=cos(angle);
      sinang=sin(angle);
      e1mcox=(1.-cosang)*u.x();
      e1mcoy=(1.-cosang)*u.y();
      e1mcoz=(1.-cosang)*u.z();
    }

    double VectorRotate::getAngle() const { return std::acos(cosang); }

    /*!
     * \brief Same as rotate(const Point) but kept for compatibility
     * \todo remove
     */
    Point VectorRotate::rotate(const Geometrybase &geo, Point p) const {
      return rotate(p);
    }

    /*!
     * Rotate point around axis specified above.
     * \param p Vector to rotate
     */
    Point VectorRotate::rotate(Point p) const {
      assert(&geoPtr!=nullptr);
      Point b(p-origin);
      geoPtr->boundary(b);           // Apply boundary conditions
      double eb=u.x()*b.x() + u.y()*b.y() + u.z()*b.z();
      p.x()=e1mcox*eb+cosang*b.x()+sinang*(u.y()*b.z() - u.z()*b.y()) + origin.x();
      p.y()=e1mcoy*eb+cosang*b.y()+sinang*(u.z()*b.x() - u.x()*b.z()) + origin.y();
      p.z()=e1mcoz*eb+cosang*b.z()+sinang*(u.x()*b.y() - u.y()*b.x()) + origin.z();
      geoPtr->boundary(p);
      return p;
    }
   
    /*!
     * \brief Quaternion rotation
     */
    void QuaternionRotate::setAxis(Geometrybase &geo, const Point &beg, const Point &end, double angle) {
      assert(&geo!=nullptr);
      assert( abs(end.len()-1)<1e-7 && "Unit vector required!!");
      geoPtr=&geo;
      double sinanghalf = sin(angle*0.5);
      q.x=end.x()*sinanghalf;
      q.y=end.y()*sinanghalf;
      q.z=end.z()*sinanghalf;
      q.w=cos(angle*0.5);
      /*generate quaternion variables for rotation*/
      double t2 =  q.w * q.x;
      double t3 =  q.w * q.y;
      double t4 =  q.w * q.z;
      double t5 = -q.x * q.x;
      double t6 =  q.x * q.y;
      double t7 =  q.x * q.z;
      double t8 = -q.y * q.y;
      double t9 =  q.y * q.z;
      double t10 = -q.z * q.z;

      d1 = t8 + t10;
      d2 = t6 - t4;
      d3 = t3 + t7;
      d4 = t4 + t6;
      d5 = t5 + t10;
      d6 = t9 - t2;
      d7 = t7 - t3;
      d8 = t2 + t9;
      d9 = t5 + t8;
    }

    Point QuaternionRotate::rotate(Point p) const {
      double newx = 2.0 * ( d1*p.x() + d2*p.y() + d3*p.z() ) + p.x();
      double newy = 2.0 * ( d4*p.x() + d5*p.y() + d6*p.z() ) + p.y();
      double newz = 2.0 * ( d7*p.x() + d8*p.y() + d9*p.z() ) + p.z();
      p.x() = newx;
      p.y() = newy;
      p.z() = newz;
      return p;
    }

    /*!
     * \param dir1 Direction of segment
     * \param halfl1 Half length of segment
     * \param Distance vector between the middle segment to point
     */
    Point mindist_segment2point(const Point &dir, double halfl, const Point &r_cm) {
      double d;
      double c = dir.dot(r_cm);
      if (c > halfl)
        d = halfl;
      else {
        if (c > -halfl) d = c;
        else d = -halfl;
      }
      return -r_cm + (dir*d);
    }

    /*!
     * \param dir1 Direction of first segment
     * \param halfl1 Half length of first segment
     * \param dir2 Direction of second segment
     * \param halfl2 Half length of second segment
     * \param r_cm Distance vector between the middle of the two segments
     */
    Point mindist_segment2segment(const Point &dir1, double halfl1,
        const Point &dir2, double halfl2, const Point &r_cm)
    {
      Point u = dir1 * (halfl1*2); //S1.P1 - S1.P0;
      Point v = dir2 * (halfl2*2); //S2.P1 - S2.P0;
      Point w = dir2*halfl2 - dir1*halfl1 - r_cm; //S1.P0-S2.P0
      double a = u.dot(u);        // always >= 0
      double b = u.dot(v);
      double c = v.dot(v);        // always >= 0
      double d = u.dot(w);
      double e = v.dot(w);
      double D = a*c - b*b;       // always >= 0
      double sc = D;
      double sN = D;
      double sD = D;      // sc = sN / sD, default sD = D >= 0
      double tc = D;
      double tN = D;
      double tD = D;      // tc = tN / tD, default tD = D >= 0

      // compute the line parameters of the two closest points
      if (D < 0.00000001) { // the lines are almost parallel
        sN = 0.0;        // force using point P0 on segment S1
        sD = 1.0;        // to prevent possible division by 0.0 later
        tN = e;
        tD = c;
      }
      else {                // get the closest points on the infinite lines
        sN = (b*e - c*d);
        tN = (a*e - b*d);
        if (sN < 0.0) {       // sc < 0 => the s=0 edge is visible
          sN = 0.0;
          tN = e;
          tD = c;
        }
        else if (sN > sD) {  // sc > 1 => the s=1 edge is visible
          sN = sD;
          tN = e + b;
          tD = c;
        }
      }

      if (tN < 0.0) {           // tc < 0 => the t=0 edge is visible
        tN = 0.0;
        // recompute sc for this edge
        if (-d < 0.0)
          sN = 0.0;
        else if (-d > a)
          sN = sD;
        else {
          sN = -d;
          sD = a;
        }
      }
      else if (tN > tD) {      // tc > 1 => the t=1 edge is visible
        tN = tD;
        // recompute sc for this edge
        if ((-d + b) < 0.0)
          sN = 0;
        else if ((-d + b) > a)
          sN = sD;
        else {
          sN = (-d + b);
          sD = a;
        }
      }
      // finally do the division to get sc and tc
      if (std::abs(sN) < 0.00000001) sc = 0.0 ;
      else sc = sN / sD;
      if (std::abs(tN) < 0.00000001) tc = 0.0 ;
      else tc = tN / tD;

      // get the difference of the two closest points
      //Vector = w + (sc * u) - (tc * v);  // = S1(sc) - S2(tc)
      //vec.x = u.x*sc + w.x - v.x*tc;
      //vec.y = u.y*sc + w.y - v.y*tc;
      //vec.z = u.z*sc + w.z - v.z*tc;
      return u*sc + w - v*tc;
    }    
      
      /*!
       \brief Calculates intersections of spherocylinder2 with a all-way patch of spherocylinder1 and return them (PSC)
       */
      int psc_intersect(const CigarParticle &part1, const CigarParticle &part2, const Point &r_cm, double intersections[5], double rcut)
      
      {
          int intrs;
          double a, b, c, d, e, x1, x2, rcut2;
          Point cm21, vec1, vec2, vec3, vec4;          
          
          intrs=0;
          rcut2=rcut*rcut;
          /*1- do intersections of spherocylinder2 with patch of spherocylinder1 at 
           cut distance C*/
          /*1a- test intersection with half planes of patch and look how far they are 
           from spherocylinder. If closer then C  we got itersection*/
          
          /* plane1 */
          /* find intersections of part2 with plane by par1 and patchsides[0] */
          intrs+=find_intersect_plane(part1,part2,r_cm,part1.patchsides[0],rcut,part1.pcanglsw,intersections);
          //	    printf("plane1 %d\n", intrs);
          /* plane2 */
          /* find intersections of part2 with plane by par1 and patchsides[1] */
          intrs+=find_intersect_plane(part1,part2,r_cm,part1.patchsides[1],rcut,part1.pcanglsw,intersections);
          
          if ( (intrs == 2 ) && (part1.pcanglsw <0) ) {
              fprintf (stderr, "ERROR: Patch is larger than 180 degrees and we are getting two segments - this hasnot been programed yet.\n\n");
              exit (1);
          }
          //	    printf("plane2 %d\n", intrs);
          
          /*1b- test intersection with cylinder - it is at distance C*/
          if (intrs < 2 )  {
              cm21=-r_cm;
              vec1=cm21.cross(part1.dir);
              vec2=part2.dir.cross(part1.dir);
              a = vec2.dot(vec2);
              b = 2*vec1.dot(vec2);
              c = -rcut*rcut + vec1.dot(vec1);
              d = b*b - 4*a*c;
              if ( d >= 0) { /*there is intersection with infinite cylinder */
                  x1 = (-b+sqrt(d))*0.5/a;/*parameter on line of SC2 determining intersection*/
                  if ((x1 >=part2.halfl) || (x1 <= -part2.halfl)) intrs+=0; /*intersection is outside sc2*/
                  else {
                      /* vectors from center os sc1 to intersection with infinite cylinder*/
                      vec1=part2.dir*x1-r_cm;
                      e = part1.dir.dot(vec1);
                      if ((e >=part1.halfl) || (e <= -part1.halfl)) intrs+=0; /*intersection is outside sc1*/
                      else {
                          intrs+=test_intrpatch(part1,vec1,part1.pcanglsw,x1,intersections);
                      }
                  }
                  if ( d > 0 ){
                      x2 = (-b-sqrt(d))*0.5/a;/*parameter on line of SC2 determining intersection*/
                      if ((x2 >=part2.halfl) || (x2 <= -part2.halfl)) intrs+=0; /*intersection is outside sc2*/
                      else {
                          vec2 = part2.dir*x2-r_cm;
                          e = part1.dir.dot(vec2);
                          if ((e >=part1.halfl) || (e <= -part1.halfl)) intrs+=0; /*intersection is outside sc1*/
                          else {
                              intrs+=test_intrpatch(part1,vec2,part1.pcanglsw,x2,intersections);
                          }
                      }
                  }
              }
          }
          //	    printf ("cylinder %d x1 %f x2 %f e %f\n", intrs, x1, x2, e);
          /*1c- test intersection with spheres at the end - it is at distace C*/
          if (intrs < 2 )  {
              /*centers of spheres*/
              /*relative to the CM of sc2*/
              vec1 =  part1.dir*part1.halfl - r_cm;
              vec2 = -part1.dir*part1.halfl - r_cm;
              
              /*sphere1*/
              a = part2.dir.dot(part2.dir);
              b = 2.0*vec1.dot(part2.dir);
              c = vec1.dot(vec1)-rcut*rcut;
              d = b*b-4*a*c;
              if (d >= 0) { /*if d<0 there are no intersections*/
                  x1= (-b + sqrt(d))*0.5/a; /*parameter on line of SC2 determining intersection*/
                  if ((x1 >=part2.halfl) || (x1 <= -part2.halfl)) intrs+=0; /*intersection is outside sc2*/
                  else {
                      vec3 = part2.dir*x1-r_cm;
                      e = part1.dir.dot(vec3);
                      if ((e >= part1.halfl) || (e <= -part1.halfl)) { /*if not intersection is inside sc1*/
                          intrs+=test_intrpatch(part1,vec3,part1.pcanglsw,x1,intersections);
                      }
                  }
                  if ( d > 0) {
                      x2= (-b - sqrt(d))*0.5/a; /*parameter on line of SC2 determining intersection*/
                      if ((x2 >=part2.halfl) || (x2 <= -part2.halfl)) intrs+=0; /*intersection is outside sc2*/
                      else {
                          vec4 = part2.dir*x2 - r_cm;
                          e = part1.dir.dot(vec4);
                          if ((e >=part1.halfl) || (e <= -part1.halfl)) { /*if not intersection is inside sc1*/
                              intrs+=test_intrpatch(part1,vec4,part1.pcanglsw,x2,intersections);
                          }
                      }
                  }
              }
              //		printf ("sphere1 %d x1 %f x2 %f e %f\n", intrs, x1, x2, e);
              /*sphere2*/
              a = part2.dir.dot(part2.dir);
              b = 2.0*vec2.dot(part2.dir);
              c = vec2.dot(vec2)-rcut*rcut;
              d = b*b-4*a*c;
              if (d >= 0) { /*if d<0 there are no intersections*/
                  x1= (-b + sqrt(d))*0.5/a; /*parameter on line of SC2 determining intersection*/
                  if ((x1 >=part2.halfl) || (x1 <= -part2.halfl)) intrs+=0; /*intersection is outside sc2*/
                  else {
                      vec3 = part2.dir*x1 - r_cm;
                      e = part1.dir.dot(vec3);
                      if ((e >=part1.halfl) || (e <= -part1.halfl)) { /*if not intersection is inside sc1*/
                          intrs+=test_intrpatch(part1,vec3,part1.pcanglsw,x1,intersections);
                      }
                  }
                  if ( d > 0 ) {
                      x2= (-b - sqrt(d))*0.5/a; /*parameter on line of SC2 determining intersection*/
                      if ((x2 >=part2.halfl) || (x2 <= -part2.halfl)) intrs+=0; /*intersection is outside sc2*/
                      else {
                          vec4 = part2.dir*x2 - r_cm;
                          e = part1.dir.dot(vec4);
                          if ((e >=part1.halfl) || (e <= -part1.halfl)) { /*if not intersection is inside sc1*/
                              intrs+=test_intrpatch(part1,vec4,part1.pcanglsw,x2,intersections);
                          }
                      }
                  }
              }
              //		printf ("sphere2 %d\n", intrs);
          }
          
          /*1d- if there is only one itersection shperocylinder ends within patch wedge 
           set as second intersection end inside patch*/
          if (intrs < 2 )  {
              /*whole spherocylinder is in or all out if intrs ==0*/
              vec1 = part2.dir*part2.halfl - r_cm;
              /*vector from CM of sc1 to end of sc2*/
              /*check is is inside sc1*/
              a=vec1.dot(part1.dir);
              vec3 = vec1 - part1.dir*a;
              b=vec3.dot(vec3);
              d = fabs(a)-part1.halfl;
              if ( d <= 0) 
                  c = b; /*is inside cylindrical part*/
              else 
                  c = d*d + b; /*is inside caps*/
              /*c is distance squared from line or end to test if is inside sc*/
              if (c < rcut2) 
                  intrs+=test_intrpatch(part1,vec1,part1.pcanglsw,part2.halfl,intersections);
              if (intrs < 2 ) {
                  vec2 = -part2.dir*part2.halfl - r_cm;
                  /*check is is inside sc1*/
                  a=vec2.dot(part1.dir);
                  vec4 = vec2 - part1.dir*a;
                  b=vec4.dot(vec4);
                  d = fabs(a) -part1.halfl;
                  if (d <= 0) 
                      c = b; /*is inside cylindrical part*/
                  else 
                      c = d*d + b; /*is inside caps*/
                  /*c is distance squared from line or end to test if is inside sc*/
                  if (c < rcut2) 
                      intrs+=test_intrpatch(part1,vec2,part1.pcanglsw,-1.0*part2.halfl,intersections);
              }
              //		    printf ("ends %d\n", intrs);
          }
          
          
          return intrs;
      }
      
      /*!
       \brief Finds if vector vector "vec" has angular intersection with a patch of part1 
       */
      int test_intrpatch(const CigarParticle &part1, Point &vec, double cospatch, 
                         double ti, double intersections[5])
      {
          double a;
          int i, intrs;
          
          intrs=0;
          /*test if we have intersection*/
          /* do projection to patch plane*/
          vec=vec_perpproject(vec,part1.dir);
          vec.normalize();
          /* test angle distance from patch*/
          a = part1.patchdir.dot(vec);
          if (a >= cospatch) {
              intrs=1;
              i=0;
              while (intersections[i] !=0) {
                  if (ti == intersections[i]) 
                      intrs=0; /* found intersection we already have -it is at boundary*/
                  i++;
              }
              if (intrs > 0) 
                  intersections[i]=ti;
          }
          
          return intrs;
      }
      
      /*!
       \brief Finds intersections of spherocylinder and plane defined by vector "w_vec" and if they are in all-way patch then returns number of them (PSC)
       */
      int find_intersect_plane(const CigarParticle &part1, const CigarParticle &part2,  
                               const Point &r_cm, const Point &w_vec, double rcut, double cospatch, double intersections[5])
      {
          int i, intrs;
          double a, c, d, ti, disti;
          Point nplane, d_vec;
          
          nplane=part1.dir.cross(w_vec);
          nplane.normalize();
          a =  nplane.dot(part2.dir);
          if (a == 0.0) intrs=0; /* there is no intersection plane and sc are paralel*/
          else {
              ti = nplane.dot(r_cm)/a;
              if ((ti  > part2.halfl ) || (ti < -part2.halfl)) intrs=0; /* there is no intersection plane sc is too short*/
              else {
                  d_vec = ti * part2.dir - r_cm; /*vector from intersection point to CM*/
                  c = d_vec.dot(w_vec);
                  if ( c * cospatch < 0) intrs=0; /* the intersection in plane is on other side of patch */
                  else {
                      d = fabs(d_vec.dot(part1.dir)) - part2.halfl;
                      if (d <= 0) disti = c*c; /*is inside cylinder*/
                      else disti = d*d + c*c; /*is inside patch*/
                      if (disti > rcut*rcut) intrs=0; /* the intersection is outside sc */
                      else {
                          intrs=1;
                          i=0;
                          while (intersections[i] !=0) {
                              if (ti == intersections[i]) intrs=0; /* found intersection we already have -it is at boundary*/
                              i++;
                          }
                          if (intrs > 0) intersections[i]=ti;
                      }
                  }
              }
              
          }
          return intrs;
      }
      
      
      /*CPSC................................................................................*/
      
      /*!
       \brief Calculates intersections of sherocylinder2 with a cylindrical patch of spherocylinder1 and return them (CPSC)
       */
      int cpsc_intersect(const CigarParticle &part1, const CigarParticle &part2,
                        const Point &r_cm, double intersections[5], double rcut)
      
      {
          int intrs;
          double a, b, c, d, e, x1, x2, rcut2;
          Point cm21, vec1, vec2, vec3, vec4;
                 
          intrs=0;
          rcut2=rcut*rcut;
          /*1- do intersections of spherocylinder2 with patch of spherocylinder1 at 
           cut distance C*/
          /*1a- test intersection with half planes of patch and look how far they are 
           from spherocylinder. If closer then C  we got itersection*/
          
          /* plane1 */
          /* find intersections of part2 with plane by par1 and part1.patchsides[0] */
          intrs+=find_intersect_planec(part1,part2,r_cm,part1.patchsides[0],rcut,part1.pcanglsw,intersections);
          //	    printf("plane1 %d\n", intrs);
          /* plane2 */
          /* find intersections of part2 with plane by par1 and part1.patchsides[1] */
          intrs+=find_intersect_planec(part1,part2,r_cm,part1.patchsides[1],rcut,part1.pcanglsw,intersections);
          
          if ( (intrs == 2 ) && (part1.pcanglsw < 0) ) {
              fprintf (stderr, "ERROR: Patch is larger than 180 degrees and we are getting two segments - this hasnot been programed yet.\n\n");
              exit (1);
          }
          //	    printf("plane2 %d\n", intrs);
          
          /*1b- test intersection with cylinder - it is at distance C*/
          if (intrs < 2 )  {
              cm21=-r_cm;
              vec1=cm21.cross(part1.dir);
              vec2=part2.dir.cross(part1.dir);
              a = vec2.dot(vec2);
              b = 2*vec1.dot(vec2);
              c = -rcut*rcut + vec1.dot(vec1);
              d = b*b - 4*a*c;
              if ( d >= 0) { /*there is intersection with infinite cylinder */
                  x1 = (-b+sqrt(d))*0.5/a; /*parameter on line of SC2 determining intersection*/
                  if ((x1 >=part2.halfl) || (x1 <= -part2.halfl)) intrs+=0; /*intersection is outside sc2*/
                  else {
                      /* vectors from center os sc1 to intersection with infinite cylinder*/
                      vec1=part2.dir*x1-r_cm;
                      e = part1.dir.dot(vec1);
                      if ((e >=part1.halfl) || (e <= -part1.halfl)) intrs+=0; /*intersection is outside sc1*/
                      else {
                          intrs+=test_intrpatch(part1,vec1,part1.pcanglsw,x1,intersections);
                      }
                  }
                  if ( d > 0 ){
                      x2 = (-b-sqrt(d))*0.5/a; /*parameter on line of SC2 determining intersection*/
                      if ((x2 >=part2.halfl) || (x2 <= -part2.halfl)) intrs+=0; /*intersection is outside sc2*/
                      else {
                          vec2 = part2.dir*x2-r_cm;
                          e = part1.dir.dot(vec2);
                          if ((e >=part1.halfl) || (e <= -part1.halfl)) intrs+=0; /*intersection is outside sc1*/
                          else {
                              intrs+=test_intrpatch(part1,vec2,part1.pcanglsw,x2,intersections);
                          }
                      }
                  }
              }
          }
          //	    printf ("cylinder %d x1 %f x2 %f e %f\n", intrs, x1, x2, e);
          /*1c- test intersection with plates at the end - it is at distace C and in wedge*/
          if (intrs < 2 )  {
              a =  part1.dir.dot( part2.dir);
              if (a == 0.0) intrs=0; /* there is no intersection plane and sc are paralel*/
              else {
                  /*plane cap1*/
                  vec1= r_cm + part1.halfl*part1.dir;
                  x1 = part1.dir.dot(vec1)/a; /*parameter on line of SC2 determining intersection*/
                  if ((x1 > part2.halfl ) || (x1 < -part2.halfl)) intrs+=0; /* there is no intersection plane sc is too short*/
                  else {
                      vec2 = x1*part2.dir - vec1; /*vector from ENDPOINT to intersection point */
                      b = vec2.dot( vec2);
                      if (b > rcut*rcut) intrs+=0;   /* the intersection is outside sc */
                      else {
                          intrs+=test_intrpatch(part1,vec2,part1.pcanglsw,x1,intersections);
                      }
                  }
                  //		    printf ("plane cap1 %d %f\n", intrs, x1);
                  /*plane cap2*/
                  vec1= r_cm - part1.halfl*part1.dir;
                  x2 = part1.dir.dot(vec1)/a; /*parameter on line of SC2 determining intersection*/
                  if ((x2  > part2.halfl ) || (x2 < -part2.halfl)) intrs+=0; /* there is no intersection plane sc is too short*/
                  else {
                      vec2 = x2*part2.dir - vec1; /*vector from ENDPOINT to intersection point */
                      b = vec2.dot( vec2);
                      if (b > rcut*rcut) intrs+=0;   /* the intersection is outside sc */
                      else {
                          intrs+=test_intrpatch(part1,vec2,part1.pcanglsw,x2,intersections);
                      }
                  }
                  //		    printf ("plane cap2 %d %f\n", intrs,x2);
                  
              }
          }
          
          /*1d- if there is only one itersection shperocylinder ends within patch wedge 
           set as second intersection end inside patch*/
          if (intrs < 2 )  {
              /*whole spherocylinder is in or all out if intrs ==0*/
              vec1 = part2.dir*part2.halfl - r_cm;
              /*vector from CM of sc1 to end of sc2*/
              /*check is is inside sc1*/
              a=vec1.dot(part1.dir);
              vec3 = vec1 - part1.dir*a;
              b=vec3.dot(vec3);
              d = fabs(a)-part1.halfl;
              if ( d <= 0) { /*is in cylindrical part*/
                  /*c is distance squared from line or end to test if is inside sc*/
                  if (b < rcut2) intrs+=test_intrpatch(part1,vec1,part1.pcanglsw,part2.halfl,intersections);
              }
              if (intrs < 2 ) {
                  vec2 = -part2.dir*part2.halfl - r_cm;
                  /*check is is inside sc1*/
                  a=vec2.dot(part1.dir);
                  vec4 = vec2 - part1.dir*a;
                  b=vec4.dot(vec4);
                  d = fabs(a) -part1.halfl;
                  if (d <= 0) {
                      /*c is distance squared from line or end to test if is inside sc*/
                      if (b < rcut2) intrs+=test_intrpatch(part1,vec2,part1.pcanglsw,-1.0*part2.halfl,intersections);
                  }
              }
              //		    printf ("ends %d\n", intrs);
          }
          
          
          return intrs;
      }
      
      /*CPSC................................................................................*/
      
      /*!
       \brief Finds intersections of plane defined by vector "w_vec" and if they are in cylindrical patch then returns number of them (CPSC)
       */
      int find_intersect_planec(const CigarParticle &part1, const CigarParticle &part2, 
                                const Point &r_cm, const Point &w_vec, double rcut, double cospatch, double intersections[5])
      {
          int i, intrs=0;
          double a, c, d, ti, disti;
          Point nplane, d_vec;
          
          nplane=part1.dir.cross(w_vec);
          nplane.normalize();
          a =  nplane.dot( part2.dir);
          if (a == 0.0) intrs=0; /* there is no intersection plane and sc are paralel*/
          else {
              ti = nplane.dot(r_cm)/a;
              if ((ti  > part2.halfl ) || (ti < -part2.halfl)) intrs=0; /* there is no intersection plane sc is too short*/
              else {
                  d_vec = ti*part2.dir - r_cm; /*vector from intersection point to CM*/
                  c = d_vec.dot( w_vec);
                  if ( c *cospatch < 0) intrs=0; /* the intersection in plane is on other side of patch */
                  else {
                      d = fabs(d_vec.dot( part1.dir)) - part2.halfl;
                      if (d <= 0) {
                          disti= c*c; /*is inside cylinder*/
                          if (disti > rcut*rcut) intrs=0; /* the intersection is outside sc */
                          else {
                              intrs=1;
                              i=0;
                              while (intersections[i] !=0) {
                                  if (ti == intersections[i]) intrs=0; /* found intersection we already have -it is at boundary*/
                                  i++;
                              }
                              if (intrs > 0) intersections[i]=ti;
                          }
                      }
                  }
              }
          }
          return intrs;
      }
      
  
      
      
      

  }//namespace geometry

}//namespace
