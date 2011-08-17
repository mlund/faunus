#include "faunus/container.h"
#include "faunus/inputfile.h"
#include "faunus/species.h"
#include "faunus/physconst.h"
#include "faunus/group.h"

namespace Faunus {
  
  double space::charge() const {
    double z=0;
    for (unsigned short i=0; i<p.size(); i++)
      z+=p[i].charge;
    return z;
  }

  bool space::check_vector() {
    bool rc=false;
    if (p.size()==trial.size())
      for (unsigned short i=0; i<p.size(); i++) {
        if (p[i].x!=trial[i].x ||
            p[i].y!=trial[i].y ||
            p[i].z!=trial[i].z ||
            p[i].charge!=trial[i].charge)
        {
          rc=false;
          break;
        } else rc=true;
      }
    if (rc==false)
      std::cerr << "# Fatal error: Particle vectors corrupted!!\n";
    return rc;
  }

  /*!
   * \param a Particle to insers
   * \param i Position in particle vector
   */
  bool space::insert(particle a, unsigned int i) {
    if (i>p.size())
      return false;
    p.insert(p.begin()+i, a);
    trial.insert(trial.begin()+i, a);

    for (int j=0; j<g.size(); j++) {  // move and expand groups if appropriate
      if ( i<g[j]->beg )
        g[j]->beg++;
      if ( i<=g[j]->end )
        g[j]->end++;
    }
    return true;
  }

  bool space::remove(unsigned int i) {
    if (i>=p.size())
      return false;
    p.erase( p.begin()+i );
    trial.erase( trial.begin()+i );
    for (int j=0; j<g.size(); j++) {  // move and reduce groups if appropriate
      if ( i<g[j]->beg )
        g[j]->beg--;
      if (i<=g[j]->end)
        g[j]->end--;
    }
    return true;
  }
  
  bool space::saveToDisk(string file) {
    if (!fout)
      fout.open( file.c_str() );
    if (fout) {
      fout.precision(10);
      fout << p.size() << endl;
      for (int i=0; i<p.size(); i++)
        fout << p[i] << endl;
      for (int i=0; i<g.size(); i++)
        fout << g[i]->beg << " " << g[i]->end << endl;
      fout.close();
      return true;
    }
    return false;
  }
  
  /*!
   * \param file Filename
   * \param resize True if the current container should be resized to match file content (default: false)
   */
  bool space::loadFromDisk(string file, bool resize) {
    unsigned int n;
    std::ifstream f(file.c_str() );
    if (f) {
      f >> n;
      if (resize==true)
        p.resize(n);
      if (n==p.size()) {
        for (int i=0; i<n; i++)
          p[i] << f;
        trial=p;
        f.close();
        std::cout << "# Read " << n << " space from " << file << endl;
        return true;
      }
      f.close();
    }
    std::cerr << "# Container data NOT read from file " << file << endl;
    return false;
  }

  double container::dist(const point &p1, const point &p2) {
    return sqrt( sqdist(p1, p2) );
  }

  point container::vdist(const point&a, const point&b) {
    return a-b;
  }

  void container::scale(point &a, const double &s) const {
  }

  string container::info() {
    double z=charge();
    std::ostringstream o;
    o << endl
      << "# SIMULATION CONTAINER:" << endl
      << "#   Number of space  = " << p.size() << endl
      << "#   Volume (AA^3)        = " << volume << endl
      << "#   Electroneutrality    = " 
      << ((abs(z)>1e-7) ? "NO!" : "Yes") << " "  << z << endl;
    return o.str();
  }

  void container::setvolume(double vol) {
    volume=vol;
  }

  bool container::saveToDisk(string file) {
    if (!fout)
      fout.open( file.c_str() );
    if (fout) {
      fout.precision(10);
      fout << p.size() << " " << getvolume() << endl;
      return space::saveToDisk(file);
    }
    return false;
  }

  /*!
   * \param file Filename
   * \param resize True if the current container should be resized to match file content (default: false)
   */
  bool container::loadFromDisk(string file, bool resize) {
    unsigned int n;
    double v;
    std::ifstream f(file.c_str() );
    if (f) {
      f >> n >> v;
      setvolume(v);
      if (resize==true)
        p.resize(n);
      if (n==p.size()) {
        for (int i=0; i<n; i++)
          p[i] << f;
        trial=p;
        f.close();
        std::cout << "# Read " << n << " space from " << file << endl;
        return true;
      }
      f.close();
    }
    std::cerr << "# Container data NOT read from file " << file << endl;
    return false;
  }

  /*!
   * In addition to the normal collision() check that checks for
   * collision with the container boundaries, it is also possible
   * to check for certain "collisions" inside the container. This
   * is done with slicecollision(). A typical usage is to restrict
   * space or molecules to certain volumes within a simulation
   * container (window sampling).
   */
  bool container::collision_internal(const particle &p) {
    return false;
  }

  //
  //--- cell container ---
  //

  void cell::setvolume(double vol) {
    volume=vol;
    setradius( pow( 3*vol/(4*pyc.pi), 1/3.) );
  }

  cell::cell(double radius) {
    setradius(radius);
  }

  cell::cell(inputfile &in)  {
    atom.load(in);
    setradius(in.getflt("cellradius"));
  }

  void cell::setradius(double radius) {
    r = radius; 
    r2 = r*r; 
    diameter = 2*r; 
    volume = (4./3.)*pyc.pi*r*r*r;
  }

  string cell::info() {
    std::ostringstream o;
    o << container::info() 
      << "#   Shape                = Spherical" << endl
      << "#   Radius               = " << r << endl;
    return o.str();
  }

  void cell::randompos(point &m) {
    double l=r2+1;
    while (l>r2) {
      m.x = slp.random_half()*diameter;
      m.y = slp.random_half()*diameter;
      m.z = slp.random_half()*diameter;
      l=m.x*m.x+m.y*m.y+m.z*m.z;
    }
  }

  void cell::boundary(point &m) const {}

  bool cell::collision(const particle &a) {
    double x,y,z;
    x=std::abs(a.x);//+a.radius;
    y=std::abs(a.y);//+a.radius;
    z=std::abs(a.z);//+a.radius;
    return ( x*x+y*y+z*z > r2 ) ? true:false;
  }

  //
  //--- cuboid container ---
  //

  cuboid::cuboid(inputfile &in) {
    atom.load(in);
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

  string cuboid::info() {
    std::ostringstream o;
    o << container::info() 
      << "#   Sidelength           = " << len.x << "x" << len.y << "x" << len.z << endl
      << "#   Slice position[x y z]= " << len_half.x-slice_max.x << "-" << len_half.x-slice_min.x << " " 
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

  bool cuboid::collision_internal(const particle &a) {
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
    return false;
  }

  //
  //--- cuboid slit container ---
  //
  
  cuboidslit::cuboidslit(inputfile &in) : cuboid(in) {
  }

  string cuboidslit::info() {
    std::ostringstream o;
    o << cuboid::info() 
      << "#   Periodicity          = slit: xy periodicity only" << endl;
    return o.str();
  }


  /*!
   * \param length Length of the cylinder
   * \param radius Radius of the cylinder
   */
  cylinder::cylinder(double length, double radius) {
    len=length;
    r=radius;
    r2=r*r;
    diameter=r*2;
    volume=2*r2*pyc.pi*len;
    halflen=len/2;
  }

  cylinder::cylinder(inputfile &in) {
    atom.load(in);
    len=in.getflt("cylinder_len", 0);
    r=in.getflt("cylinder_radius", 0);
    r2=r*r;
    diameter=r*2;
    volume=2*r2*pyc.pi*len;
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

  bool cylinder::collision(const particle &a) {
    return 
      ( a.x*a.x+a.y*a.y>r2 || ( a.z<-halflen || a.z>halflen ) ) ? true:false;
  }

  string cylinder::info() {
    std::ostringstream o;
    o << container::info()
      << "#   Shape                = Cylindrical" << endl
      << "#   Lenght               = " << len <<endl
      << "#   Radius               = " << r << endl;
    return o.str();
  }

#ifdef HYPERSPHERE
  const double hypersphere::pi=3.141592654;

  hypersphere::hypersphere(inputfile &in) : cell(in) {
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
