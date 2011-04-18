#include <faunus/group.h>
#include "faunus/container.h"
#include "faunus/species.h"
#include "faunus/iobabel.h"
#include "faunus/physconst.h"
#include "faunus/io.h"

namespace Faunus {

  //
  // G R O U P
  //

  group::group(int i) {
    set(-1,-1);     //is empty
    title="GROUP";
    name="(noname)";//is nameless
  }

  void group::set(int first, int last) {
    beg=first;
    end=last;
  }

  int group::size() const { return (beg<0 || end<0) ? 0 : end-beg+1; }

  /*! \brief Calculates total charge
   *  \return Charge number. The charge is also stored in
   *          cm.charge.
   */
  double group::charge(const vector<particle> &p) {
    double z=0;
    for (int i=beg; i<=end; i++)
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

  int group::random() { return beg + slp.rand() % size(); }

  string group::info() {
    std::ostringstream o;
    o << endl
      << "# " << title << ": " << name << endl
      << "#   Number of particles      = " << size() << " [" << beg << "," << end << "]" << endl
      << "#   Mass center              = " << cm << endl;
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
    for (int i=beg; i<=end; i++) {
      cm += p[i] * p[i].mw;
      sum += p[i].mw;
    }
    cm=cm*(1./sum);
    cm_trial = cm;
    assert(sum>0);
    assert(beg<=end);
    return cm;
  }

  /*!
   * This mass-center routine obeys periodic boundaries (if only)
   * by first centering the first particle in (0,0,0), calc. CM
   * and move it back.
   *
   * \author Mikael Lund
   * \date Dec. 2007, Prague
   * \todo Use masscenter(con,p) function implemented below.
   */
  point group::masscenter(const container &con) {
    double sum=0;
    cm.clear();
    point t, o = con.p.at(beg); // set origo to first particle
    for (int i=beg; i<=end; i++) {
      t = con.p[i]-o;        // translate to origo
      con.boundary(t);       // periodic boundary (if any)
      cm += t * con.p[i].mw;
      sum += con.p[i].mw; 
    }
    cm=cm*(1./sum) + o;
    con.boundary(cm);
    cm_trial = cm;
    assert(sum>0);
    return cm;
  }

  /*!
   * This mass-center routine obeys periodic boundaries (if any)
   * by first centering the first particle in (0,0,0), calc. CM
   * and moving it back. Note that this function does NOT update the
   * cm and cm_trial vectors!
   *
   * \author Chris Evers
   * \date Nov. 2010, Lund
   */
  point group::masscenter(const container &con, const vector<particle> &p) const {
    double sum=0;
    point tcm;               // temporary mass center
    point t, o = p.at(beg);  // set origo to first particle
    for (int i=beg; i<=end; i++) {
      t = p[i]-o;            // translate to origo
      con.boundary(t);       // periodic boundary (if any)
      tcm += t * p[i].mw;
      sum += p[i].mw; 
    }
    tcm=tcm*(1./sum) + o;
    con.boundary(tcm);
    assert(sum>0);
    return tcm;
  }

  void group::undo(particles &par) {
    cm_trial = cm;
    for (int i=beg; i<=end; i++) {
#ifndef HYPERSPHERE
      par.trial[i].x = par.p[i].x;
      par.trial[i].y = par.p[i].y;
      par.trial[i].z = par.p[i].z;
#else
      par.trial[i].z1 = par.p[i].z1;
      par.trial[i].z2 = par.p[i].z2;
      par.trial[i].z3 = par.p[i].z3;
      par.trial[i].z4 = par.p[i].z4;
#endif
    }
  }

  void group::accept(particles &par) {
    cm = cm_trial;
    int i,ilen=end+1;
    for (i=beg; i<ilen; i++) {
#ifndef HYPERSPHERE
      par.p[i].x = par.trial[i].x;
      par.p[i].y = par.trial[i].y;
      par.p[i].z = par.trial[i].z;
#else
      par.p[i].z1 = par.trial[i].z1;
      par.p[i].z2 = par.trial[i].z2;
      par.p[i].z3 = par.trial[i].z3;
      par.p[i].z4 = par.trial[i].z4;
#endif
    }
  }

  void group::add(container &con, vector<particle> v, bool collision) {
    beg=con.p.size();
    for (int i=0; i<v.size(); i++) {
      con.boundary(v[i]);
      con.push_back(v[i]);
      end=beg+i;
    }
    masscenter(con);                // calc. mass center

    // test for overlap w. other particles
    if (collision==true) {
      move(con, -cm);               // translate to origo (0,0,0)
      accept(con); 
      point a;
      while (overlap(con)==true) {  // No overlap allowed
        con.randompos(a);           // find a random, vacent position
        a=a-cm;
        con.boundary(a);
        move(con,a);
        accept(con);
      }
    }
  }

  void group::add(container &con, unsigned char id, int n) {
    particle a=atom(id);
    if (beg<0)
      beg=con.p.size();
    for (int i=0; i<n; i++) {
      do con.randompos(a);
      while (con.overlap(a)==true); 
      end = con.push_back(a)-1;
    }
  }

  int group::displace(container &c, point dp) {
    int i=random();
#ifndef HYPERSPHERE
    c.trial[i].x = c.p[i].x + dp.x*slp.random_half();
    c.trial[i].y = c.p[i].y + dp.y*slp.random_half();
    c.trial[i].z = c.p[i].z + dp.z*slp.random_half();
    c.boundary(c.trial[i]);
#else
    double dangle=dp.x;
    double nfi=2.*acos(-1.)*slp.random_one();
    double nrho=sqrt((slp.random_one()-1.)*sin(dangle)*sin(dangle)+1.);
    double nomega=(2.*slp.random_one()-1.)*acos(cos(dangle)/nrho);
    c.trial[i].move(nrho,nomega,nfi);
#endif
    return i;
  }

  /*!
   * \param par Container class
   * \param c ...by adding this vector to all particles
   */
  void group::move(container &par, point c) {
    for (int i=beg; i<=end; i++) {
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
    int i,j;
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

  int group::count(vector<particle> &p, unsigned char id) {
    int i,n=0;
    for (i=0; i<p.size(); i++)
      if (p[i].id==id)
        n++;
    return n;
  }

  int group::nummolecules() { return size(); }

  unsigned short group::numhydrophobic( vector<particle> &p ) {
    int n=0;
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

  void group::invert(vector<particle> &p, point &iv) {
    for (int i=beg; i<=end; i++) {
      p[i].x=-p[i].x+2*iv.x;
      p[i].y=-p[i].y+2*iv.y;
      p[i].z=-p[i].z+2*iv.z;
    }
  }

  bool group::swap(container &c, int pos) {
    if (beg!=pos && size()>0) {
      std::swap_ranges(c.p.begin()+beg, c.p.begin()+end, c.p.begin()+pos);
      std::swap_ranges(c.trial.begin()+beg, c.trial.begin()+end, c.trial.begin()+pos);
      end=pos+size()-1;
      beg=pos;
      return true;
    }
    return false;
  }
  
  /*! Save all charges in group to disk
   */
  bool group::saveCharges(string filename, vector<particle> &p){
    io fio;
    std::ostringstream o;
    o.precision(12);
    o << size() << endl;
    for (int i=beg; i<=end; i++) {
      string name=atom[p[i].id].name;
      o << i << "   " << name << "   " << p[i].charge << endl;
    }
    return fio.writefile(filename, o.str());
  }

  /*! Load charges in group from disk
   */
  bool group::loadCharges(string filename, vector<particle> &p){
    std::ifstream f(filename.c_str());
    int n,i;
    double q;
    string s;
    if (f) {
      f >> n;
      for (int k=0; k<n; k++) {
        f >> i >> s >> q;
        p[i].charge=q;
      }
      f.close();
      cout << "# Loaded polymer charges from " << filename << endl;
      return true;
    } else
      std::cerr << "!! charge file not loaded !!" << endl;
    return false;
  }

  /*****************************
   *
   *          S A L T
   *
   *****************************/
  salt::salt(unsigned char cat, unsigned char an)
  {
    title="MOBILE SALT PARTICLES";
    cation=cat;
    anion=cat;
    name="MOBILE-SALT";
  }

  string salt::info(container &con) {
    std::ostringstream o;
    float c=1./pyc.Nav/1e-27;
    int nan=count(con.p, anion),
        ncat=count(con.p, cation);
    o << group::info();
    o << "#   Cation (id,z,r,N,conc) = "
      << atom[cation].name << ", " << atom[cation].charge << ", "
      << atom[cation].radius << ", "
      << ncat << ", " << ncat/con.getvolume()*c << endl
      << "#   Anion (id,z,r,N,conc)  = "
      << atom[anion].name << ", " << atom[anion].charge << ", "
      << atom[anion].radius << ", "
      << nan << ", " << nan/con.getvolume()*c << endl;
    return o.str();
  }

  void salt::isobaricmove(container &con, double newlen) {
    for (int i=beg; i<=end; i++)
      con.scale(con.trial[i], newlen);
  }

  /*!
   * This searches the inputfile object for the
   * keywords "nion#" and "tion#" and if found
   * tries to insert the salt particles at random
   * positions.
   */
  void salt::add(container &con, inputfile &in) {
    int n=1, npart;
    unsigned char id;
    while (n<15) {
      std::ostringstream nion, tion;
      nion << "nion" << n;
      tion << "tion" << n++;
      npart = in.getint(nion.str(), 0);
      id = atom( in.getstr( tion.str() ) ).id;
      if (npart!=0)
        group::add(con, id, npart ); // add particles
      if (atom[id].charge>0)
        cation=id;
      if (atom[id].charge<0)
        anion=id;
    }
  }

  //--------------- MACROMOLECULE ---------------
  macromolecule::macromolecule() { title="MACROMOLECULE"; }

  string macromolecule::info() {
    char w=10;
    std::ostringstream o;
    o << group::info();
    o << std::left;
    if (Q.cnt>0)
      o << "#   Mean charge              = "  << setw(w) << Q.avg() << endl
        << "#    Fluctuat. <Z2>-<Z>2     = "  << setw(w) << Q2.avg()-pow(Q.avg(),2) << endl;
    if (rg2.cnt>0)
      o << "#   Rms radius of gyration   = " << setw(w) << sqrt(rg2.avg()) << endl
        << "#    Fluctuat. <Rg2>-<Rg>2   = " << setw(w) << rg2.avg()-pow(rg.avg(),2) << endl;
    if (ree2.cnt>0)
      o << "#   Rms end-to-end distance  = " << setw(w) << sqrt(ree2.avg()) << endl
        << "#    Fluctuat. <Ree2>-<Ree>2 = " << setw(w) << ree2.avg()-pow(ree.avg(),2) << endl;
    if (dip.cnt>0)
      o << "#     <mu> <mu2>-<mu>2       = " << setw(w) << dip.avg() << setw(w) << dip2.avg()-pow(dip.avg(),2) << endl
        << "#     Dipole vector          = " << setw(w) << mu << endl;
      //o << "#    Radius                = " << cm.radius << endl;
    return o.str();
  }

  string macromolecule::info(container &con) {
    if (cm.radius<=0)
      vradius(con.p);
    std::ostringstream o;
    o << info();
    o << "#   Current charge         = " << charge(con.p) << endl
      << "#   Hydrophobic particles  = " << numhydrophobic(con.p) << endl
      << "#   Radius                 = " << cm.radius;
    return o.str();
  }

  int macromolecule::nummolecules() { return 1; }

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
    for (int i=beg; i<=end; ++i) {
      mu.x += (p[i].x-cm.x) * p[i].charge;
      mu.y += (p[i].y-cm.y) * p[i].charge;
      mu.z += (p[i].z-cm.z) * p[i].charge;
    }
    a=mu.len();
    dip+=a;    // avg. scalar
    dip2+=a*a; // avg. scalar squared
    return a;
  }

  double macromolecule::dipole(const container &con) {
    masscenter(con);
    particle tmp;
    mu.clear();
    for (int i=beg; i<=end; ++i) {
      tmp=con.p[i];
      tmp.x = tmp.x - cm.x;
      tmp.y = tmp.y - cm.y;
      tmp.z = tmp.z - cm.z;
      con.boundary(tmp);
      mu += tmp*tmp.charge;
    }
    double a=mu.len();
    dip+=a;    // avg. scalar
    dip2+=a*a; // avg. scalar squared
    return a;
  }

  double macromolecule::radius(vector<particle> &p) {
    double r,max=0;
    masscenter(p);
    for (int i=beg; i<=end; i++) {
      r=p[i].dist(cm) + p[i].radius;
      if (r>max)
        max=r;
    }
    cm.radius = max;
    return max;
  } 

  double macromolecule::gradius(vector<particle> &p) {
    double r=0;
    int cnt=0;
    masscenter(p);
    for (int i=beg; i<=end; i++) {
      r+=p[i].dist(cm) + p[i].radius;
      cnt++;
    }
    r/=float(cnt);
    cm.radius = r;
    return r;
  }
  
  /*!
   * This squared radius of gyration routine is mass weighted and obeys periodic boundaries 
   * (if any) by following the previous defined masscenter method
   */    
  double macromolecule::sqmassgradius(container &c) {
    double x=0;
    double r2=0;
    double sum=0;
    point t, o;
//    masscenter(c);                      // recalculate center of mass
    for (int i=beg; i<=end; i++) {
      t = c.p[i]-cm;                    // vector to center of mass
      c.boundary(t);                    // periodic boundary (if any)
      r2 = c.sqdist(t,o);               // squared distance to cm
      x += r2 * c.p[i].mw;              // m*r^2
      sum += c.p[i].mw;                 // total mass
    }
    x=x*(1./sum);
    rg2+=x;                             // add to average squared radius of gyration
    rg+=sqrt(x);                        // add to average radius of gyration
    return x;                           // obtained squared radius of gyration
  }

  /*!
   * This squared radius of gyration routine does not update any averages,
   * is mass weighted and obeys periodic boundaries (if any) 
   * by following the previous defined masscenter method
   */
  double macromolecule::sqmassgradius(container &c, vector<particle> &p) {
    double x=0;
    double r2=0;
    double sum=0;
    point t, o;
    for (int i=beg; i<=end; i++) {
      t = p[i]-cm;                      // vector to center of mass
      c.boundary(t);                    // periodic boundary (if any)
      r2 = c.sqdist(t,o);               // squared distance to cm
      x += r2 * c.p[i].mw;              // m*r^2
      sum += c.p[i].mw;                 // total mass
    }
    x=x*(1./sum);
    return x;                           // obtained squared radius of gyration
  }

  double macromolecule::sqend2enddistance(container &c) {
    double x=0;
    x=c.sqdist(c.p[beg], c.p[end]);     // squared distance between begin and end
    ree2+=x;                            // add to average squared end-to-end distance
    ree+=sqrt(x);                       // add to average end-to-end distance
    return x;                           // obtained squared radius of gyration
  }

  double macromolecule::vradius(vector<particle> &p) {
    double v=0,r;
    for (int i=beg; i<=end; i++) {
      v+=1.333333*3.14159265*p[i].radius*p[i].radius*p[i].radius;
    }
    r=pow( 0.75 / 3.14159265 * v , 1./3.);
    cm.radius = r;
    return r;
  }

  double macromolecule::charge(const vector<particle> &p) {
    double z=group::charge(p);
    Q+=z;
    Q2+=z*z;
    return z;
  }

  double macromolecule::getcharge(const vector<particle> &p) {
    return group::charge(p);
  }

  void macromolecule::center(container &con) {
    move(con, -cm);
    accept(con);
  }

  void macromolecule::zmove(container &par, double dz) {
    cm_trial.z = cm.z + dz;
    par.boundary(cm_trial);
    for (int i=beg; i<=end; i++) {
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
   * \param par Container
   * \param u Vector to rotate around
   * \param angle ..by this many degrees (rad)
   * \param dr Displacement parameter
   * \todo This could be done more elegant using quaternion parameters (Frenkel p49)
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
    for (int i=beg; i<=end; i++) {
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
    //cm_trial = masscenter(par.trial); // calling masscenter() is no good in trial move as it sets cm_trial=cm
    cm_trial = cm + q;
    par.boundary(cm_trial);
  }

  /*!
   * \param par Container
   * \param cr Center of rotation, center of mass of seed
   * \param u Vector to rotate around
   * \param angle ..by this many degrees (rad)
   * \warning Ej bijektiv avbildning med periodiska randvilkor. Anvand bindningsvektorer
   */
  void macromolecule::rotate(container &par, point cr, point u, double angle) {
    point b,c;
    double cosang, sinang, eb;
    double e1mcox, e1mcoy, e1mcoz;

    cosang=cos(angle);
    sinang=sin(angle);
    e1mcox=(1.-cosang)*u.x;
    e1mcoy=(1.-cosang)*u.y;
    e1mcoz=(1.-cosang)*u.z;
//    c=par.vdist(cm,cr);
//    for (int i=beg; i<=end; i++) {        // Unwrap the boundary conditions 
  //    par.trial[i]=cm+par.vdist(par.p[i],cm);
    //}

    for (int i=beg; i<=end; i++) {
      b.x=par.trial[i].x - cr.x;              // translate to cr...
      b.y=par.trial[i].y - cr.y;
      b.z=par.trial[i].z - cr.z;
      par.boundary(b);                    // Apply boundary conditions
      // No, wait a second
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
   * Simultaneous rotation and translation 
   *
   * \param par Container
   * \param dr Displacement parameter
   * \param angle ..by this many degrees (rad)
   * \warning Ej bijektiv avbildning med periodiska randvilkor. Anvand bindningsvektorer
   */
  void macromolecule::transrot(container &par, double dr, double angle) {
    // Let us first translate...
    point c;
    c.ranunit(slp); 
    c=c*dr; //random displacement vector
    for (int i=beg; i<=end; i++) {
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
    for (int i=beg; i<=end; i++) {
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

  //void macromolecule::add(container &con, inputfile &in ) { }

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

  int molecules::random() {
    return (group::random()-beg)/numatom;
  }

  group molecules::operator[](int i) {
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
          if (c.clash(c.p[j], p[i])==true) {
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

  polymer::polymer() {}

#ifdef BABEL
  bool polymer::babeladd(container &c, inputfile &in) {
    name=in.getstr("polymer");
    return babeladd(c, in, name);
  }

  bool polymer::babeladd(container &c, inputfile &in, int &num) {
    string parameter;
    std::stringstream buffer;
    buffer << "polymer" << num+1;
    parameter = buffer.str();
    name=in.getstr(parameter,"none");
    if ( name == "none" )
      name=in.getstr("polymer");
    return babeladd(c, in, name);
  }

  bool polymer::babeladd(container &c, inputfile &in, string &name) {
    nb.clear();
    iobabel ob;
    ob.read(name);
    if (ob.p.size()>0) {
      add(c,ob.p);
      move(c,-cm);
      accept(c);
      for (int i=beg; i<=end; i++) {
        nb.push_back( ob.neighbors(i-beg) );
        for (int j=0; j<nb[i-beg].size(); j++)
          nb[i-beg][j]+=beg;
      }
    double d2 = .75 * pow(c.getvolume(),2./3.);  // Squared half box space diagonal
    for (int k=0; k<c.p.size()-1; k++) 
      for (int l=1; l<c.p.size(); l++) {
        double r2 = pow(c.p[k].x-c.p[l].x,2)+pow(c.p[k].y-c.p[l].y,2)+pow(c.p[k].z-c.p[l].z,2); 
        if ( r2>d2 ) {
          cout << "! WARNING! Box to small for polymer! " << endl;
          cout << "! Half box space diagnonal   = " << sqrt(d2) << endl;
          cout << "! Distance between monomers  = " << sqrt(r2) << endl;
          return false;
        }
      }
      return true;
    }
    return false;
  }
#endif

  vector<int> polymer::neighbors(int i) const {
    return nb.at(i-beg);
  }

  bool polymer::addbond(int i, int j) {
    if (i>=beg && i<=end)
      if (j>=beg && j<=end)
        if (areneighbors(i,j)==false) {
          if ( (i-beg)<=nb.size() ) nb.resize(i-beg+1);
          if ( (j-beg)<=nb.size() ) nb.resize(j-beg+1);
          nb[i-beg].push_back(j);
          nb[j-beg].push_back(i);
          return true;
        }
    return false;
  }

  //!< is j a neighbor to polymer atom i?
  bool polymer::areneighbors(int i, int j) const {
    //    return ( std::find( nb[i-beg].begin(), nb[i-beg].end(), j) != nb[i-end].end() ) ? true : false;
    //    std::cout <<"particel "<<i<<" with "<<nb[i-beg].size()<<"n, looking for "<<j<<std::endl;
    for (int k=0; k<nb[i-beg].size(); k++)
      if (nb[i-beg][k]==j)
        return true;
    return false;
  }

  string polymer::info() {
    std::ostringstream o;
    o << macromolecule::info();
    if (nb.size()<100.) {                   // Only print connectivity if the number of bounds < 100
      o << "#   Connectivity:" << std::endl;
      for (int i=0; i<nb.size(); i++) {
        o << "#     Atom " << i+beg << ": ";
        for (int j=0; j<nb[i].size(); j++)
          o << nb[i][j] << " ";
        o << std::endl;
      }
    }
    return o.str();
  }

  /*!
   * This will generate a set of VMD commands that will
   * create bonds between neighboring monomers. This can
   * conveniently be put into a tcl script file. VMD can
   * subsequently be started with:\n
   *
   * vmd confout.pqr -e vmdbonds.tcl
   *
   */
  string polymer::getVMDBondScript() {
    std::ostringstream o;
    for (int i=0; i<nb.size(); i++) {
      o << "set sel [atomselect top \"index " << i+beg << "\"]" << endl
        << "$sel setbonds [list { ";
      for (int j=0; j<nb[i].size(); j++) {
        o << nb[i][j] << " ";
      }
      o << "}]" << endl;
    }
    return o.str(); 
  }

  popscmembrane::popscmembrane() {
    name.assign("Mixed Phospholipid Membrane (POPS/C)");
  }

  string popscmembrane::info() {
    int ns=pops.size(), nc=popc.size(), ntot=ns+nc;
    std::ostringstream o;
    o << group::info()
      << "#   Area per headgroup     = " << headarea<<" A^2" << endl
      << "#   POPS:POPC ratio        = " << ns*100./ntot << ":" << nc*100./ntot << endl
      << "#   Number of lipids       = " << ntot
      << " (" << ns << " pops / " << nc << " popc)" << endl;
    return o.str();
  }

  string popscmembrane::info(slit &c) {
    double z=charge(c.p);
    std::ostringstream o;
    o.precision(4);
    o << popscmembrane::info()
      << "#   Total charge           = " << z << endl;
    if (z!=0)
      o << "#   Charge density         = " << c.xyarea/z << " A^2/e" << endl;
    return o.str();
  }

  void popscmembrane::load(inputfile &in, slit &con) {
    //Some objects and variables
    int bbeg=con.p.size();
    int npoly, npolys, npolyc, graftp;
    double gstep;
    particle a;
    a.mw=1.;
    vector<particle> Pops, Popc;
    polymer p;
    Pops.clear(), Popc.clear();
    popc.clear(), pops.clear();
    //Prepare particle vector of popc/s
    a.radius=2.5, /*a.id='PLG',*/ a.charge=0, a.z=-con.zlen*0.5;
    Pops.push_back(a), Popc.push_back(a);
    /*a.id="PL1",*/ a.charge=-1, a.z=-con.zlen*0.5+5;
    Pops.push_back(a), Popc.push_back(a);
    /*a.id="PL2",*/ a.charge=1,  a.z=-con.zlen*0.5+10;
    Pops.push_back(a), Popc.push_back(a);
    /*a.id="PL1",*/ a.charge=-1, a.z=-con.zlen*0.5+15;
    Pops.push_back(a);
    //Prepare the membrane and grid
    scratio =in.getflt("scratio",-1);
    headarea=in.getflt("headarea", -1);
    npoly   =int(con.xyarea/headarea);
    npolys  =int(double(npoly)*scratio);
    npolyc  =npoly-npolys;
    graftp  =int(sqrt(double(npoly))+1.);     //No. gridpoints along a axis
    gstep   =sqrt(con.xyarea)/(double(graftp));
    double x,y;
    x=y=-con.len_half;
    int cnt=0;
    vector<int> nb;
    for (int i=0; i<graftp; i++) {
      Pops[0].x=Pops[1].x=Pops[2].x=Pops[3].x=x+i*gstep;
      Popc[0].x=Popc[1].x=Popc[2].x=x+i*gstep;
      for (int j=0; j<graftp; j++) {
        //        nb.clear();
        //        p.nb.clear();
        ++cnt;
        Pops[0].y=Pops[1].y=Pops[2].y=Pops[3].y=y+j*gstep;
        Popc[0].y=Popc[1].y=Popc[2].y=y+j*gstep;
        if(cnt<=npolys) {  //Push back and set neigbour lists
          p.add(con, Pops);
          p.nb.clear();
          nb.clear();
          p.nb.resize(4);
          nb.push_back(1+p.beg);
          p.nb[0]=nb;
          nb.clear();
          nb.push_back(p.beg),   nb.push_back(p.beg+2);
          p.nb[1]=nb;
          nb.clear();
          nb.push_back(1+p.beg), nb.push_back(3+p.beg);
          p.nb[2]=nb;
          nb.clear();
          nb.push_back(2+p.beg);
          p.nb[3]=nb;
          pops.push_back(p);
        }
        if(cnt>npolys && cnt<=npoly) {
          p.add(con, Popc);
          p.nb.clear();
          nb.clear();
          p.nb.resize(3);
          nb.push_back(1+p.beg);
          p.nb[0]=nb;
          nb.clear();
          nb.push_back(p.beg),   nb.push_back(p.beg+2);
          p.nb[1]=nb; 
          nb.clear();
          nb.push_back(1+p.beg);
          p.nb[2]=nb;
          popc.push_back(p);
        }
      }
    }
    //    for (int i=0; i<popc.size(); i++)
    //      for (int j=0; j<popc[i].nb.size(); j++)
    //        for (int k=0; j<popc[i].nb[j].size(); k++)
    //          std::cout <<popc[i].nb[j][k]<<std::endl;
    beg=bbeg;
    end=con.p.size()-1;
  }

  string popscmembrane::getVMDBondScript() {
    std::ostringstream o;
    for (int i=0; i<pops.size(); i++)
      o << pops[i].getVMDBondScript();
    for (int i=0; i<popc.size(); i++)
      o << popc[i].getVMDBondScript();
    return o.str();
  }

#ifdef BABEL
  glu3::glu3(container &con, inputfile &in) {
    name="GLU3 DENDRIMER";
    chains.babeladd(con, in);
    beg=chains.beg, end=chains.end;
    masscenter(con);
    core.beg=beg, core.end=beg+51;
    core.name="PORPHYRIN CORE";
    chains.beg=beg+52;
    chains.nb.erase(chains.nb.begin(), chains.nb.begin()+52);
    chains.name="GLUTAMIC CHAINS";
    /*    cout <<chains.nb.size() <<" size of nb, "<<chains.size()<<" size of chains"<<endl
          <<"glu3 starts and ends at "<<beg<<" "<<end<<" and chains at "<<chains.beg<<" "<<chains.end<<endl
          <<"porphoryn core starts and ends at "<<core.beg<<" "<<core.end<<endl;
          */
  }
  string glu3::info() {
    std::ostringstream o;
    o << macromolecule::info() <<core.info()
      << chains.info(); 
    return o.str();
  }
#endif

#ifdef HYPERSPHERE
  /*!
   * \param c Hypersphere container
   * \param v Vector of particles to add
   * \param collision True if collision check is needed. (Default = False).
   */
  void hypermolecule::add(container &c, vector<particle> v, bool collision) {
  }

  /*!
   * \param c Container
   * \param p point to add to all coordinates in group
   */
  void hypermolecule::move(container &c, point p) {
  }

  void hypermolecule::rotate(container &c, double drot, double dp) {
  }
#endif
}//namespace
