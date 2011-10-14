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
  point Group::masscenter(const Space &con) const {
    double sum=0;
    point cm,t,o = con.p.at(beg); // set origo to first particle
    for (int i=beg; i<=end; ++i) {
      t = con.p[i]-o;        // translate to origo
      con.geo->boundary(t);       // periodic boundary (if any)
      cm += t * con.p[i].mw;
      sum += con.p[i].mw; 
    }
    cm=cm*(1/sum) + o;
    con.geo->boundary(cm);
    assert(sum>0);
    return cm;
  }
  
  point Group::dipolemoment(const Space &con) const {
    point d;
    return d;
  }
  
  void Group::rotate(Space &con, const point &v, double) {
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
  void Group::translate(Space &con, const point &p) {
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

  Atomic::Atomic() {
    id=ATOMIC;
  }

  Atomic::Atomic(Space &spc, InputMap &in) {
    add(spc,in);
  }

  void Atomic::add(Space &spc, InputMap &in) {
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

  Atomic& Atomic::operator<<(std::istream &in) {
    Group::operator<<(in);
    return *this;
  }

  void Atomic::scale(Space&, double) {
  }


  /* -----------------------------------*
   *             MOLECULAR
   * -----------------------------------*/
  Molecular::Molecular() {
    id=MOLECULAR;
  }

  std::ostream& Molecular::write(std::ostream &o) const {
    Group::write(o);
    o << " " << cm;
    return o;
  }

  Molecular& Molecular::operator<<(std::istream &in) {
    Group::operator<<(in);
    cm.operator<<(in);
    return *this;
  }

  void Molecular::translate(Space &con, const point &p) {
    Group::translate(con, p);
    cm_trial=cm+p;
    con.geo->boundary(cm_trial);
  }

  void Molecular::accept(Space &s) {
    Group::accept(s);
    cm=cm_trial;
  }

  void Molecular::scale(Space&, double) {
  }


}//namespace
