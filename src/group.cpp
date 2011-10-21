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

  Group::~Group() {}

  Group::Group(int first, int last) : beg(first), end(last) {
    id=GROUP;
  }

  int Group::size() const {
    return (beg<0 || end<0) ? 0 : end-beg+1;
  }

  bool Group::find(int i) const {
    return (i>=beg && i<=end) ? true : false;
  }

  /*! \brief Calculates total charge
   *  \return Charge number. The charge is also stored in
   *          cm.charge.
   */
  double Group::charge(const p_vec &p) const {
    double z=0;
    for (int i=beg; i<=end; ++i)
      z+=p[i].charge;
    return z;
  }

  Group& Group::operator+=(const Group& g) {
    if (g.beg==-1 && g.end==-1)   // if added Group is empty
      return (*this);
    if (beg==-1 && end==-1)       // if this is empty
      beg=g.beg;
    end=g.end;
    return (*this);
  }

  const Group Group::operator+(const Group& g) const {
    // TODO: this should probably be done using operator += like this:
    // return Group(*this) += g;
    // (move extra functionality to +=)
    Group o;
    if (beg<g.beg) {
      o.beg=beg;
      o.end=g.end;
    } else {
      o.beg=g.beg;
      o.end=end;
    };
    o.name = this->name + " + " + g.name;
    assert( o.size()==this->size()+g.size() ); // debug
    if (o.size()!=this->size()+g.size())
      std::cout << "# Warning: Added Groups are not continous!\n";
    return o;
  }

  bool Group::operator==(const Group& g) const {
    return (*this == g);
  }
  
  std::ostream& Group::write(std::ostream &o) const {
    o << beg << " " << end;
    return o;
  }

  std::ostream& operator<<(std::ostream &o, Group &g) {
    return g.write(o);
  }
  
  Group& Group::operator<<(std::istream &in) {
    in >> beg >> end;
    return *this;
  }

  string Group::info() {
    std::ostringstream o;
    return o.str();
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
  Point Group::_massCenter(const Space &spc) const {
    double sum=0;
    Point cm,t,o = spc.p.at(beg); // set origo to first particle
    for (int i=beg; i<=end; ++i) {
      t = spc.p[i]-o;        // translate to origo
      spc.geo->boundary(t);       // periodic boundary (if any)
      cm += t * spc.p[i].mw;
      sum += spc.p[i].mw;
    }
    assert(sum>0); // Masses must be assigned!
    cm=cm*(1/sum) + o;
    spc.geo->boundary(cm);
    return cm;
  }

  Point Group::massCenter(const Space &spc) const {
    return _massCenter(spc);
  }
  
  Point Group::dipolemoment(const Space &con) const {
    Point d;
    return d;
  }
  
  void Group::rotate(Space &con, const Point &v, double) {
  }
  
  void Group::scale(Space &con, double) {
  }

  void Group::undo(Space &s) {
    for (int i=beg; i<=end; ++i)
      s.trial[i]=s.p[i];
  }

  void Group::accept(Space &s) {
    for (int i=beg; i<=end; ++i)
      s.p[i] = s.trial[i];
  }

  /*!
   * \param par Container class
   * \param c ...by adding this vector to all particles
   */
  void Group::translate(Space &con, const Point &p) {
    for (int i=beg; i<=end; ++i) {
      con.trial[i] = con.p[i] + p;
      con.geo->boundary( con.trial[i] );
    }
  }

  int Group::random() const {
    return beg + slp_global.rand() % size();
  }

  /* -----------------------------------*
   *             ATOMIC
   * -----------------------------------*/

  GroupAtomic::GroupAtomic() {
    id=ATOMIC;
  }

  GroupAtomic::GroupAtomic(Space &spc, InputMap &in) {
    add(spc,in);
  }

  void GroupAtomic::add(Space &spc, InputMap &in) {
    beg=spc.p.size();
    end=beg-1;
    int n=1, npart;
    do {
      std::ostringstream nion("nion"), tion("tion"), dpion("dpion");
      nion << "nion" << n;
      tion << "tion" << n;
      dpion<< "dpion"<< n++;
      npart = in.get(nion.str(), 0);
      if (npart>0) {
        short id = atom[ in.get<string>(tion.str(), "UNK") ].id;
        atom[id].dp = in.get(dpion.str(), 0.);
        spc.insert( atom[id].name, npart);
        end+=npart;
      } else break;
    } while (npart>0);
    if (beg>end)
      beg=end=-1;
    else
      spc.enroll(*this);
  }

  GroupAtomic& GroupAtomic::operator<<(std::istream &in) {
    Group::operator<<(in);
    return *this;
  }

  void GroupAtomic::scale(Space&, double) {
  }


  /* -----------------------------------*
   *             MOLECULAR
   * -----------------------------------*/
  GroupMolecular::GroupMolecular() {
    id=MOLECULAR;
  }

  std::ostream& GroupMolecular::write(std::ostream &o) const {
    Group::write(o);
    o << " " << cm;
    return o;
  }

  GroupMolecular& GroupMolecular::operator<<(std::istream &in) {
    Group::operator<<(in);
    cm.operator<<(in);
    return *this;
  }
  
  void GroupMolecular::rotate(Space &spc, const Point &endpoint, double angle) {
    assert( spc.geo->dist(cm, cm_trial)<1e-9 );   // debug. Is trial mass center in sync?
    vrot.setAxis(*spc.geo, cm, endpoint, angle);  // rotate around line between mass center and point
    for (int i=beg; i<=end; i++)
      spc.trial[i] = vrot.rotate( *spc.geo, spc.p[i] ); // boundary conditions already taken care of
  }
  
  void GroupMolecular::translate(Space &spc, const Point &p) {
    assert( spc.geo->dist(cm,cm_trial)<1e-9 );   // debug. Is trial mass center in sync?
    Group::translate(spc, p);
    cm_trial=cm+p;
    spc.geo->boundary(cm_trial);
  }

  void GroupMolecular::accept(Space &s) {
    Group::accept(s);
    cm=cm_trial;
  }

  void GroupMolecular::scale(Space&, double) {
  }

}//namespace
