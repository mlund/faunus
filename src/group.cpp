#include "faunus/group.h"

namespace Faunus {

  //---------------- GROUP -----------------------
  group::group(int i) {
    set(-1,-1);     //is empty
    title="GROUP";
    name="(noname)";//is nameless
  }
  void group::set(short int first, short int last) {
    beg=first;
    end=last;
  }
  short int group::size() const { return (beg<0 || end<0) ? 0 : end-beg+1; }

  /*! \brief Calculates total charge
   *  \return Charge number. The charge is also stored in
   *          cm.charge.
   */
  double group::charge(const vector<particle> &p) {
    double z=0;
    for (short i=beg; i<=end; i++)
      z+=p[i].charge;
    cm.charge=z;
    return z;
  }

  group& group::operator+=(const group& g) {
    if (g.beg==-1 && g.end==-1)   // if added group is empty
      return (*this);
    if (beg==-1 && end==-1)       // if this is empty
      beg=g.beg;
    end=g.end;
    return (*this);
  }

  const group group::operator+(const group& g) const {
    // TODO: this should probably be done using operator += like this:
    // return group(*this) += g;
    // (move extra functionality to +=)
    group o;
    if (beg<g.beg) {
      o.beg=beg;
      o.end=g.end;
    } else {
      o.beg=g.beg;
      o.end=end;
    };
    o.name = this->name + " + " + g.name;
    if (o.size()!=this->size()+g.size())
      std::cout << "# Warning: Added groups are not continous!\n";
    return o;
  }

  bool group::operator==(const group& g) const {
    // TODO: perhaps there could even be a real comparison (?)
    return (*this == g);
  }

  short int group::random() { return beg + rand() % size(); }

  string group::info() {
    std::ostringstream o;
    o << endl
      << "# " << title << ": " << name << endl
      << "#   Range                  = [" << beg << "," << end << "]" << " " << size() << endl
      << "#   Mass center            = " << cm << endl;
    return o.str();
  }

  bool group::find(unsigned int i) const { return (i>=beg && i<=end) ? true : false; }

  /*!
   * Recalcute mass center of the particles
   * in a group. The result is stored in the
   * group CM placeholder as well as returned
   * by the function.
   * \warning Doesn't consider periodic boundaries!
   * \note Use masscenter(container &) instead.
   */
  point group::masscenter(const vector<particle> &p) {
    cm.clear();
    double sum=0;
    for (short i=beg; i<=end; i++) {
      cm += p[i] * p[i].mw;
      sum += p[i].mw;
    }
    cm=cm*(1./sum);
    cm_trial = cm;
    return cm;
  }

  /*!
   * This mass-center routine obeys periodic boundaries (if only)
   * by first centering the first particle in (0,0,0), calc. CM
   * and move it back.
   *
   * \author Mikael Lund
   * \date Dec. 2007, Prague
   */
  point group::masscenter(container &con) {
    double sum=0;
    cm.clear();
    point t, o = con.p[beg]; // set origo to first particle
    for (unsigned short i=beg; i<=end; i++) {
      t = con.p[i]-o;        // translate to origo
      con.boundary(t);       // periodic boundary (if any)
      cm += t * con.p[i].mw;
      sum += con.p[i].mw; 
    }
    cm=cm*(1./sum) + o;
    con.boundary(cm);
    cm_trial = cm;
    return cm;
  }


  void group::undo(particles &par) {
    cm_trial = cm;
    for (short i=beg; i<=end; i++) {
      par.trial[i].x = par.p[i].x;
      par.trial[i].y = par.p[i].y;
      par.trial[i].z = par.p[i].z;
    }
  }
  void group::accept(particles &par) {
    cm = cm_trial;
    unsigned short i,ilen=end+1;
    for (i=beg; i<ilen; i++) {
      par.p[i].x = par.trial[i].x;
      par.p[i].y = par.trial[i].y;
      par.p[i].z = par.trial[i].z;
    }
  }
  void group::add(container &par, vector<particle> v, bool collision) {
    beg=par.p.size();
    for (unsigned short i=0; i<v.size(); i++) {
      par.boundary(v[i]);
      par.push_back(v[i]);
      end=beg+i;
    }
    masscenter(par);                // calc. mass center

    // test for overlap w. other particles
    if (collision==true) {
      move(par, -cm);               // translate to origo (0,0,0)
      accept(par); 
      point a;
      while (overlap(par)==true) {  // No overlap allowed
        par.randompos(a);           // find a random, vacent 
        move(par, a);               // position.
        accept(par);
      }
    }
  }
  void group::add(container &con, particle::type id, short n) {
    particle a=con.atom(id);
    if (beg<0)
      beg=con.p.size();
    for (unsigned short i=0; i<n; i++) {
      do con.randompos(a);
      while (con.overlap(a)==true); 
      end = con.push_back(a)-1;
    }
  }

  unsigned short group::displace(container &c, double dp) {
    unsigned short i=random();
    c.trial[i].x = c.p[i].x + dp*slp.random_half();
    c.trial[i].y = c.p[i].y + dp*slp.random_half();
    c.trial[i].z = c.p[i].z + dp*slp.random_half();
    c.boundary(c.trial[i]);
    return i;
  }

  /*!
   * \param par Container class
   * \param c ...by adding this vector to all particles
   * \param k Keyword to specify if the move should be automatically accepted (default no).
   */
  void group::move(container &par, point c) {
    for (unsigned short i=beg; i<=end; i++) {
      par.trial[i].x = par.p[i].x + c.x;
      par.trial[i].y = par.p[i].y + c.y;
      par.trial[i].z = par.p[i].z + c.z;
      par.boundary(par.trial[i]);
    }
    cm_trial.x = cm.x + c.x;
    cm_trial.y = cm.y + c.y;
    cm_trial.z = cm.z + c.z;
    par.boundary(cm_trial);
  }

  bool group::overlap(container &c) {
    double s;
    unsigned short i,j;
    for (i=beg; i<=end; i++) {
      if (c.collision(c.p[i])==true)  // check cell collision
        return true;
      for (j=0; j<beg; j++) {
        s=c.p[i].radius+c.p[j].radius;
        if ( c.sqdist(c.p[i],c.p[j]) < s*s )
          return true;
      }
      for (j=end+1; j<c.p.size(); j++) {
        s=c.p[i].radius+c.p[j].radius;
        if ( c.sqdist(c.p[i],c.p[j]) < s*s )
          return true;
      }
    }
    return false;
  }
  short int group::count(vector<particle> &p, particle::type id) {
    short int i,n=0;
    for (i=0; i<p.size(); i++)
      if (p[i].id==id)
        n++;
    return n;
  }
  unsigned short group::nummolecules() { return size(); }
  unsigned short group::numhydrophobic( vector<particle> &p ) {
    unsigned short n=0;
    for (int i=beg; i<=end; i++)
      if (p[i].hydrophobic==true) n++;
    return n;
  }

  bool group::swap(container &c, group &g) {
    int n=size();
    if (n!=g.size() ||
        find(g.beg) ||
        find(g.end) ) return false;
    for (int i=0; i<n; i++) {
      std::swap( c.p[ beg+i ], c.p[ g.beg+i ] );
      std::swap( c.trial[ beg+i ], c.trial[ g.beg+i ] );
    }
    std::swap(beg, g.beg);
    std::swap(end, g.end);
    return true;
  }

  /*****************************
   *
   *          S A L T
   *
   *****************************/
  salt::salt(particle::type cat, particle::type an)
  {
    title="MOBILE SALT PARTICLES";
    cation=cat;
    anion=cat;
  }
  string salt::info(container &con) {
    std::ostringstream o;
    float c=1./6.022e23/1e-27;
    short nan=count(con.p, anion),
          ncat=count(con.p, cation);
    o << group::info();
    o << "#   Cation (id,z,r,N,conc) = "
      << con.atom[cation].name << ", " << con.atom[cation].charge << ", "
      << con.atom[cation].radius << ", "
      << ncat << ", " << ncat/con.getvolume()*c << endl
      << "#   Anion (id,z,r,N,conc)  = "
      << con.atom[anion].name << ", " << con.atom[anion].charge << ", "
      << con.atom[anion].radius << ", "
      << nan << ", " << nan/con.getvolume()*c << endl;
    return o.str();
  }
  void salt::isobaricmove(container &con, double newlen) {
    for (unsigned short i=beg; i<=end; i++)
      con.scale(con.trial[i], newlen);
  }

  /*!
   * This searches the inputfile object for the
   * keywords "nion#" and "tion#" and if found
   * tries to insert the salt particles at random
   * positions.
   */
  void salt::add(container &con, inputfile &in) {
    short n=1, npart;
    particle::type id;
    while (n<4) {
      std::ostringstream nion, tion;
      nion << "nion" << n;
      tion << "tion" << n++;
      npart = in.getint(nion.str(), 0);
      id = con.atom( in.getstr( tion.str() ) ).id;
      if (npart!=0)
        group::add(con, id, npart ); // add particles
      if (con.atom[id].charge>0)
        cation=id;
      if (con.atom[id].charge<0)
        anion=id;
    }
  }

  //--------------- MACROMOLECULE ---------------
  macromolecule::macromolecule() { title="MACROMOLECULE"; }
  string macromolecule::info() {
    std::ostringstream o;
    o << group::info();
    if (Q.cnt>0)
      o << "#  Charge fluctuations:" << endl
        << "#    <Z> <Z^2>-<Z>^2       = " << Q.avg()<<" "<<Q2.avg()-pow(Q.avg(),2)<<endl;
    if (dip.cnt>0)
      o << "#    <mu> <mu^2>-<mu>^2    = " << dip.avg()<<" "<<dip2.avg()-pow(dip.avg(),2)<<endl
        << "#    Dipole vector         = " << mu << endl;
    return o.str();
  }
  string macromolecule::info(container &con) {
    std::ostringstream o;
    o << info();
    o << "#   Current charge         = " << charge(con.p) << endl
      << "#   Hydrophobic particles  = " << numhydrophobic(con.p) << endl
      << "#   Radius                 = " << radius(con.p) << endl;
    return o.str();
  }
  unsigned short macromolecule::nummolecules() { return 1; }
  macromolecule& macromolecule::operator=(const group& g) {
    beg=g.beg;
    end=g.end;
    name=g.name;
    cm=g.cm;
    cm_trial=g.cm_trial;
    return (*this);
  }
  /*!
   * Recalulates the dipole moment of a
   * group and stores the result in the
   * group dipole moment placeholder.
   * Origo is (0,0,0).
   *
   * \todo Implement PBC!
   */
  double macromolecule::dipole(vector<particle> &p)
  {
    mu.clear();
    double a;
    for (unsigned short i=beg; i<=end; ++i) {
      mu.x += (p[i].x-cm.x) * p[i].charge;
      mu.y += (p[i].y-cm.y) * p[i].charge;
      mu.z += (p[i].z-cm.z) * p[i].charge;
    }
    a=mu.len();
    dip+=a;    // avg. scalar
    dip2+=a*a; // avg. scalar squared
    return a;
  }
  double macromolecule::radius(vector<particle> &p) {
    double r,max=0;
    masscenter(p);
    for (short i=beg; i<=end; i++) {
      r=p[i].dist(cm) + p[i].radius;
      if (r>max)
        max=r;
    }
    cm.radius = max;
    return max;
  } 
  double macromolecule::charge(const vector<particle> &p) {
    double z=group::charge(p);
    Q+=z;
    Q2+=z*z;
    return z;
  }
  void macromolecule::center(container &con) {
    move(con, -cm);
    accept(con);
  };
  void macromolecule::zmove(container &par, double dz) {
    cm_trial.z = cm.z + dz;
    par.boundary(cm_trial);
    for (short int i=beg; i<=end; i++) {
      par.trial[i].z = par.p[i].z + dz;
      par.boundary(par.trial[i]);
    };
  }
  void macromolecule::rotate(container &par, double drot, double dp) {
    point u;
    double r=2;
    while (r > 1.) { //Generate a random unit vector
      u.x=slp.random_one();
      u.y=slp.random_one();
      u.z=slp.random_one();
      r=sqrt(u.x*u.x+u.y*u.y+u.z*u.z);
    }
    u.x=u.x/r;
    u.y=u.y/r; 
    u.z=u.z/r;
    rotate(par, u, drot*slp.random_half(), dp);
  }

  /*!
   * \param g Group to rotate
   * \param u Vector to rotate around
   * \param angle ..by this many degrees (rad)
   * \param k Keyword to specify automatic acceptance
   * \toto This could be done more elegant using quaternion parameters (Frenkel p49)
   */
  void macromolecule::rotate(container &par, point u, double angle, double dr) {
    point b, q;
    q.x=dr*slp.random_half();
    q.y=dr*slp.random_half();
    q.z=dr*slp.random_half();
    double cosang, sinang, eb;
    double e1mcox, e1mcoy, e1mcoz;
    cosang=cos(angle);
    sinang=sin(angle);
    e1mcox=(1.-cosang)*u.x;
    e1mcoy=(1.-cosang)*u.y;
    e1mcoz=(1.-cosang)*u.z;
    for (short i=beg; i<=end; i++) {
      b.x=par.p[i].x - cm.x;              // translate to origo...
      b.y=par.p[i].y - cm.y;
      b.z=par.p[i].z - cm.z;
      par.boundary(b);                    // Apply boundary conditions
      eb=u.x*b.x + u.y*b.y + u.z*b.z;
      par.trial[i].x=e1mcox*eb+cosang*b.x+sinang*(u.y*b.z - u.z*b.y) + cm.x + q.x;
      par.trial[i].y=e1mcoy*eb+cosang*b.y+sinang*(u.z*b.x - u.x*b.z) + cm.y + q.y;
      par.trial[i].z=e1mcoz*eb+cosang*b.z+sinang*(u.x*b.y - u.y*b.x) + cm.z + q.z;
      par.boundary(par.trial[i]);
      
    }
    cm_trial = cm + q;
    par.boundary(cm_trial);
  }

  /*!
   * \param g Group to rotate around arbitrary point 
   * \param cr Center of rotation, center of mass of seed
   * \param u Vector to rotate around
   * \param angle ..by this many degrees (rad)
   * \param k Keyword to specify automatic acceptance
   * \warning Ej bijektiv avbildning med periodiska randvilkor
   * \warning Anvand bindningsvektorer
   */
  void macromolecule::rotate(container &par, point cr, point u, double angle) {
    point b;
    double cosang, sinang, eb;
    double e1mcox, e1mcoy, e1mcoz;

    cosang=cos(angle);
    sinang=sin(angle);
    e1mcox=(1.-cosang)*u.x;
    e1mcoy=(1.-cosang)*u.y;
    e1mcoz=(1.-cosang)*u.z;
    for (short i=beg; i<=end; i++) {
      b.x=par.p[i].x - cr.x;              // translate to cr...
      b.y=par.p[i].y - cr.y;
      b.z=par.p[i].z - cr.z;
      par.boundary(b);                    // Apply boundary conditions
      eb=u.x*b.x + u.y*b.y + u.z*b.z;
      par.trial[i].x=e1mcox*eb+cosang*b.x+sinang*(u.y*b.z - u.z*b.y) + cr.x;
      par.trial[i].y=e1mcoy*eb+cosang*b.y+sinang*(u.z*b.x - u.x*b.z) + cr.y;
      par.trial[i].z=e1mcoz*eb+cosang*b.z+sinang*(u.x*b.y - u.y*b.x) + cr.z;
      par.boundary(par.trial[i]);
    }
    b.x=cm.x - cr.x;              // translate to cr...
    b.y=cm.y - cr.y;
    b.z=cm.z - cr.z;
    par.boundary(b);                    // Apply boundary conditions
    eb=u.x*b.x + u.y*b.y + u.z*b.z;
    cm_trial.x=e1mcox*eb+cosang*b.x+sinang*(u.y*b.z - u.z*b.y) + cr.x;
    cm_trial.y=e1mcoy*eb+cosang*b.y+sinang*(u.z*b.x - u.x*b.z) + cr.y;
    cm_trial.z=e1mcoz*eb+cosang*b.z+sinang*(u.x*b.y - u.y*b.x) + cr.z;
    par.boundary(cm_trial);
  }

  /*!
   * \param g Consecutice rotation and translation 
   * \param c displacement vector
   * \param u Vector to rotate around
   * \param angle ..by this many degrees (rad)
   * \param k Keyword to specify automatic acceptance
   * \warning Ej bijektiv avbildning med periodiska randvilkor
   * \warning Anvand bindningsvektorer
   */
  void macromolecule::transrot(container &par, double dr, double angle) {
   // Let us first translate...
    point c;
    c.ranunit(slp); 
    c=c*dr; //random displacement vector
    for (unsigned short i=beg; i<=end; i++) {
      par.trial[i].x = par.p[i].x + c.x;
      par.trial[i].y = par.p[i].y + c.y;
      par.trial[i].z = par.p[i].z + c.z;
      par.boundary(par.trial[i]);
    }
    cm_trial.x = cm.x + c.x;
    cm_trial.y = cm.y + c.y;
    cm_trial.z = cm.z + c.z;
    par.boundary(cm_trial);

    //... then we rotate
    c.ranunit(slp);
    double ang;
    ang=angle*slp.random_half();
    point b;
    double cosang, sinang, eb;
    double e1mcox, e1mcoy, e1mcoz;

    cosang=cos(ang);
    sinang=sin(ang);
    e1mcox=(1.-cosang)*c.x;
    e1mcoy=(1.-cosang)*c.y;
    e1mcoz=(1.-cosang)*c.z;
    for (short i=beg; i<=end; i++) {
      b.x=par.trial[i].x - cm_trial.x;              // center cm_trial...
      b.y=par.trial[i].y - cm_trial.y;
      b.z=par.trial[i].z - cm_trial.z;
      par.boundary(b);                    // Apply boundary conditions
      eb=c.x*b.x + c.y*b.y + c.z*b.z;
      par.trial[i].x=e1mcox*eb+cosang*b.x+sinang*(c.y*b.z - c.z*b.y) + cm_trial.x; // move back in to laboratory frame
      par.trial[i].y=e1mcoy*eb+cosang*b.y+sinang*(c.z*b.x - c.x*b.z) + cm_trial.y;
      par.trial[i].z=e1mcoz*eb+cosang*b.z+sinang*(c.x*b.y - c.y*b.x) + cm_trial.z;
      par.boundary(par.trial[i]);
    }
  }

  void macromolecule::add(container &con, inputfile &in ) { }

  /*! \warning Assumes a cubic box!
   */
  void macromolecule::isobaricmove(container &con, double newlen) {
    double oldvol=con.getvolume(); // store original volume
    con.scale(cm_trial, newlen);   // scale center of mass (cm)
    point CM_trial=cm_trial;       // copy trial cm
    move(con, -cm );               // move to origo
    con.setvolume(pow(newlen,3));  // apply trial volume
    for (int i=beg;i<=end;i++) {
      con.trial[i]+= CM_trial;     // move all particles to new cm
      con.boundary(con.trial[i]);  // respect boundary conditions
    }
    cm_trial=CM_trial;
    con.setvolume(oldvol);         // restore original volume
  }

  //--------------- MOLECULES -----------------
  molecules::molecules(unsigned short n) {
    title="MOLECULAR ARRAY";
    numatom=n;
  }
  short int molecules::random() {
    return (group::random()-beg)/numatom;
  }
  group molecules::operator[](unsigned short i) {
    sel.beg=beg+i*numatom;
    sel.end=sel.beg+numatom-1;
    return sel;
  }
  string molecules::info() {
    std::ostringstream o;
    o << macromolecule::info()
      << "#   Atoms per molecule     = " << numatom << endl 
      << "#   Solvent molecules      = " << size()/double(numatom) << endl; //check!!
    return o.str();
  }
  void molecules::add(container &c, vector<particle> &p, short step) {
    bool col;
    beg=end=c.p.size();
    for (int i=0; i<p.size(); i+=step) {
      col=c.collision(p[i]);
      if (col==false)
        for (int j=0; j<c.p.size(); j++)
          if (c.collision(c.p[j], p[i])==true) {
            col=true;
            break;
          }
      if (col==false)
        for (int k=i; k<i+step; k++) {
          c.trial.push_back(p[k]);
          end=c.trial.size()-1;
        }
    }
    c.p=c.trial;
  }
  vector<int> molecules::pick(int m) {
    vector<int> n;
    n.clear();
    //find m pointers
    n.push_back(int(size()*slp.random_one()/numatom));
    int c=0, o=0;
    while (n.size()<m) {
      c=1;
      o=int(size()*slp.random_one()/numatom);
      for (int i=0; i<n.size(); i++)
        if (n[i]==o)
          c=0;
      if (c==1)
        n.push_back(o);
    }
    //sort the pointers 
    for (int i =0; i<m-1; i++)
      for (int j=i+1; j<m; j++)
        if(n[i]>n[j]) {
          o=n[j];
          n[j]=n[i];
          n[i]=o;
        }
    if (m!=n.size())
      std::cerr << "# MOLECULES::PICK Iterator out of sync with 'm'"<<endl;
    return n;
  }

  //--------------- CHAIN -----------------
  chain::chain(container &con, int n, particle::type sort, double r) { graftpoint=false; 
    req=r;
    add(con, n, sort);
  }
  void chain::add(container &con, int n, particle::type sort) { // no check of PARTICLE overlap
    point ru;
    particle a=con.atom(sort);
    beg=con.p.size();
    con.randompos(a);
    end=con.push_back(a)-1;
    ru.ranunit(slp);
    a+=ru*req;
    for (int i=1; i<n; i++) {
      do {
        a=a-ru*req;
        ru.ranunit(slp);
        a+=ru*req;
      } while (con.overlap(a));
      end = con.push_back(a)-1;
    }
  }

  chain::chain(container &con, int n, particle::type sort, double r, point &gp) { graftpoint=true; 
    GP=&gp;
    req=r;
    addgrafted(con, n, sort, gp);
  }

  void chain::addgrafted(container &con, int n, particle::type sort, point &gp) { // no check of PARTICLE overlap
    point ru;
    particle a=con.atom(sort);
    beg=con.p.size();
    con.randompos(a);
    a=gp;
    ru.ranunit(slp);
    a+=ru*req;
    for (int i=0; i<n; i++) {
      do {
        a=a-ru*req;
        ru.ranunit(slp);
        a+=ru*req;
      } while (con.overlap(a));
      end = con.push_back(a)-1;
    }
  }

  double chain::monomerenergy(vector<particle> &p, short i) {
    double u=0;
    //the first ?
    if (i==beg) {
      u+=quadratic( p[i], p[i+1] );
      if (graftpoint)
        u+=quadratic( p[i], *GP );
      return u;
    }
    //the last ?
    if (i==end)
      return quadratic( p[i], p[i-1] );
    //otherwise...
    return quadratic(p[i], p[i+1]) + quadratic(p[i], p[i-1]);
  }
  double chain::internalenergy(vector<particle> &p) {
    double u=0;
    for (short i=beg; i<end; i++)
      u+=quadratic( p[i], p[i+1]);
    if (graftpoint)
      u+=quadratic( p[beg], *GP );
    return u;
  }

  //-------------------- PLANE -----------------
  //! By default the plane lies in the XY plane at
  //! z=0. This plane can be put in XZ and YZ by changing
  //! the plane vector. To shift the plane, modify the
  //! offset vector.
  planarsurface::planarsurface() {
    plane.x=1;
    plane.y=1;
    plane.z=0;
  }
  void planarsurface::project(point &p) { p=p*plane + offset; }

  //--------------------- PHOSPHOMEMBRANE -------------------
  double zwittermembrane::selfenergy(particles &par) {
    double u=0;
    for (short i=beg; i<end; i=i+2)
      u+=pairpot(par.p[i], par.p[i+1]);
    return u;
  }
  short zwittermembrane::mate(short i) { return (i%2==0) ? i+1 : i-1; }

  //! Even particles can move only in the plane while odd
  //! particles can move in any direction. The latter
  //! corresponds to flexible N(CH3)+ groups attached to
  //! phosphate groups. 
  unsigned short zwittermembrane::displace(container &c, double d) {
    unsigned short i=group::displace(c,d);
    if (i%2==0)
      project(c.trial[i]); // keep even particles in the plane
    return i;
  } 
  // Brush
  /*
     brush::brush
     t*/      

}
