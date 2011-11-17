#include <faunus/geometry.h>
#include <faunus/point.h>
#include <faunus/species.h>
#include <faunus/physconst.h>
#include <faunus/group.h>
#include <faunus/average.h>
#include <faunus/space.h>
#include <faunus/inputfile.h>

namespace Faunus {

  /* -----------------------------------*
   *                GROUP
   * -----------------------------------*/

  Group::Group(int front, int back) : myrange(front,back-front+1) {
    id=GROUP;
  }

  Group::~Group() {}

  /*
  int Group::size() const {
    return (beg<0 || last<0) ? 0 : last-beg+1;
  }
  
  bool Group::empty() const {
    if (beg>=0)
      if (last>=beg)
	return false;
    return true;
  }
*/
  bool Group::find(int i) const {
    return (i>=front() && i<=back()) ? true : false;
  }

  /*! \brief Calculates total charge
   *  \return Charge number. The charge is also stored in
   *          cm.charge.
   */
  double Group::charge(const p_vec &p) const {
    double z=0;
    for (auto i : *this)
      z+=p[i].charge;
    return z;
  }

  Group& Group::operator+=(const Group& g) {
    assert(!"Unimplemented!");
    /*
    if (g.beg==-1 && g.last==-1)   // if added Group is empty
      return (*this);
    if (beg==-1 && last==-1)       // if this is empty
      beg=g.beg;
    last=g.last;
    */
    return (*this);
  }

  const Group Group::operator+(const Group& g) const {
    // TODO: this should probably be done using operator += like this:
    // return Group(*this) += g;
    // (move extra functionality to +=)
    assert(!"Unimplemented!");
    Group o;
    /*
    if (beg<g.beg) {
      o.beg=beg;
      o.last=g.last;
    } else {
      o.beg=g.beg;
      o.last=last;
    };
    o.name = this->name + " + " + g.name;
    assert( o.size()==this->size()+g.size() ); // debug
    if (o.size()!=this->size()+g.size())
      std::cout << "# Warning: Added Groups are not continous!\n";
      */
    return o;
  }

  bool Group::operator==(const Group& g) const {
    return (*this == g);
  }
  
  std::ostream& Group::write(std::ostream &o) const {
    o << front() << " " << back() << " " << cm;
    return o;
  }

  std::ostream& operator<<(std::ostream &o, Group &g) {
    return g.write(o);
  }
  
  Group& Group::operator<<(std::istream &in) {
    int front;
    int back;
    in >> front >> back;
    setrange(front,back);
    assert( size()==back-front+1 && "Problem with Group range");
    cm.operator<<(in);
    return *this;
  }

  string Group::info() {
    std::ostringstream o;
    return o.str();
  }

  /*!
   * This mass-center routine obeys periodic boundaries (if any)
   * by first centering the middle particle in (0,0,0), calc. CM
   * and move it back.
   *
   * \author Mikael Lund
   * \date Dec. 2007, Prague
   * \warning Intra-molecular distances must not exceed half the box size for cubouid geometry.
   * \todo Implement assertion to catch failure when molecule is bigger than half the box size.
   */
  Point Group::_massCenter(const Space &spc) const {
    double sum=0;
    Point cm,t,o = spc.p.at(front()+(back()-front())/2);  // set origo to middle particle
    for (auto i : *this) {
      t = spc.p[i]-o;              // translate to origo
      spc.geo->boundary(t);        // periodic boundary (if any)
      cm += t * spc.p[i].mw;
      sum += spc.p[i].mw;
    }
    assert(sum>0 && "Group has zero mass. Did you forget to assign atom weights?");
    cm=cm*(1/sum) + o;
    spc.geo->boundary(cm);
    return cm;
  }

  Point Group::setMassCenter(const Space &spc) {
    cm=massCenter(spc);
    cm_trial=cm;
    return cm;
  }

  Point Group::massCenter(const Space &spc) const {
    assert( &spc!=NULL );
    assert( size()>0 );
    assert( back() < (int)spc.p.size() );
    return _massCenter(spc);
  }
  
  Point Group::dipolemoment(const Space &con) const {
    assert(!"Unimplemented!");
    Point d;
    return d;
  }
  
  void Group::rotate(Space &spc, const Point &endpoint, double angle) {
    assert(!"Unimplemented!");
  }
  
  void Group::scale(Space &s, double newvol) {
    cm_trial=cm;
    s.geo->scale(cm_trial, newvol);
    for (auto i : *this)
      s.geo->scale( s.trial[i], newvol);
  }

  void Group::undo(Space &s) {
    for (auto i : *this)
      s.trial[i]=s.p[i];
    cm_trial=cm;
  }

  void Group::accept(Space &s) {
    for (auto i : *this)
      s.p[i] = s.trial[i];
    cm=cm_trial;
  }

  /*!
   * \param spc Simulation Space
   * \param p Vector to translate with
   */
  void Group::translate(Space &spc, const Point &p) {
    for (auto i : *this) {
      spc.trial[i] += p;
      spc.geo->boundary( spc.trial[i] );
    }
    cm_trial = cm + p;
    spc.geo->boundary( cm_trial );
  }

  int Group::random() const {
    int i = front() + slp_global.rand() % size();
    assert(i>=front() && i<=back() && "Generated random element out of range!");
    return i;
  }

  /* -----------------------------------*
   *             ATOMIC
   * -----------------------------------*/

  GroupAtomic::GroupAtomic() {
    id=ATOMIC;
    property.insert(ATOMIC);
  }

  GroupAtomic::GroupAtomic(Space &spc, InputMap &in) {
    id=ATOMIC;
    property.insert(ATOMIC);
    add(spc,in);
  }

  void GroupAtomic::add(Space &spc, InputMap &in) {
    setfront( spc.p.size() );
    int size=0;
    int n=1, npart;
    do {
      std::ostringstream nion("nion"), tion("tion"), dpion("dpion"), aion("aion");
      nion << "nion" << n;
      tion << "tion" << n;
      dpion<< "dpion"<< n;
      aion << "aion" << n++; //activity
      npart = in.get(nion.str(), 0);
      if (npart>0) {
        short id = atom[ in.get<string>(tion.str(), "UNK") ].id;
        atom[id].dp = in.get(dpion.str(), 0.);
        atom[id].activity = in.get(aion.str(), 0.);
        spc.insert( atom[id].name, npart);
        size+=npart;
      } else break;
    } while (npart>0);
    if (size>0) {
      resize(size);
      spc.enroll(*this);
      setMassCenter(spc);
    }
  }

  GroupAtomic& GroupAtomic::operator<<(std::istream &in) {
    Group::operator<<(in);
    return *this;
  }

  /* -----------------------------------*
   *             MOLECULAR
   * -----------------------------------*/
  GroupMolecular::GroupMolecular() {
    id=MOLECULAR;
    property.insert(MOLECULAR);
  }

  std::ostream& GroupMolecular::write(std::ostream &o) const {
    Group::write(o);
    //o << " " << cm;
    return o;
  }

  GroupMolecular& GroupMolecular::operator<<(std::istream &in) {
    Group::operator<<(in);
    //cm.operator<<(in);
    return *this;
  }
  
  void GroupMolecular::translate(Space &spc, const Point &p) {
    assert( spc.geo->dist(cm,massCenter(spc))<1e-9 );   // debug. Is trial mass center in sync?
    Group::translate(spc, p);
    cm_trial=cm+p;
    spc.geo->boundary(cm_trial);
  }

  void GroupMolecular::rotate(Space &spc, const Point &endpoint, double angle) {
    assert( spc.geo->dist(cm,massCenter(spc) )<1e-6 );      // debug. Is mass center in sync?

    cm_trial = cm;
    vrot.setAxis(*spc.geo, cm, endpoint, angle);            // rotate around line between mass center and point
    for (auto i : *this)
      spc.trial[i] = vrot.rotate( *spc.geo, spc.trial[i] ); // boundary conditions already taken care of

    assert( spc.geo->dist(cm_trial, massCenter(spc))<1e-9 && "Rotation messed up mass center. Is the box too small?");
  }

  void GroupMolecular::scale(Space &s, double newvol) {
    assert( s.geo->dist(cm, massCenter(s))<1e-6);
    assert( s.geo->dist(cm, cm_trial)<1e-7);

    Point newcm=cm;
    s.geo->scale(newcm, newvol);      // scale cm to newcm
    translate(s,-cm);                 // move to origo

    double oldvol=s.geo->getVolume(); // store original volume
    s.geo->setVolume(newvol);         // apply trial volume

    for (auto i : *this) {
      s.trial[i] += newcm;            // move all particles to new cm
      s.geo->boundary( s.trial[i] );  // respect boundary conditions
    }
    cm_trial=newcm;
    s.geo->setVolume(oldvol);         // restore original volume
  }

}//namespace
