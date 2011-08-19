#include <faunus/point.h>
#include <faunus/container.h>
#include <faunus/species.h>
#include <faunus/physconst.h>
#include <faunus/group.h>

namespace Faunus {

  /* -----------------------------------*
   *                GROUP
   * -----------------------------------*/

  group::group() {
    beg=end=-1;
  }

  int group::len() const {
    return (beg<0 || end<0) ? 0 : end-beg+1;
  }

  /*! \brief Calculates total charge
   *  \return Charge number. The charge is also stored in
   *          cm.charge.
   */
  double group::charge(const vector<particle> &p) const {
    double z=0;
    for (int i=beg; i<=end; ++i)
      z+=p[i].charge;
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
    if (o.len()!=this->len()+g.len())
      std::cout << "# Warning: Added groups are not continous!\n";
    return o;
  }

  bool group::operator==(const group& g) const { return (*this == g); }
  
  std::ostream& group::write(std::ostream &o) const {
    o << beg << " " << end;
    return o;
  }

  std::ostream& operator<<(std::ostream &o, group &g) { return g.write(o); }
  
  group & group::operator<<(std::istream &in) {
    in >> beg >> end;
    return *this;
  }


  string group::info() {
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
  point group::masscenter(const container &con) const {
    double sum=0;
    point cm,t,o = con.p.at(beg); // set origo to first particle
    for (int i=beg; i<=end; i++) {
      t = con.p[i]-o;        // translate to origo
      con.boundary(t);       // periodic boundary (if any)
      cm += t * con.p[i].mw;
      sum += con.p[i].mw; 
    }
    cm=cm*(1./sum) + o;
    con.boundary(cm);
    assert(sum>0);
    return cm;
  }
  
  point group::dipolemoment(const container &con) const {
    point d;
    return d;
  }
  
  void group::rotate(container &con, const point &v, double) {
  }
  
  void group::scale(container &con, double) {
  }

  void group::undo(space &s) {
    for (int i=beg; i<=end; i++) {
      s.trial[i].x = s.p[i].x;
      s.trial[i].y = s.p[i].y;
      s.trial[i].z = s.p[i].z;
    }
  }

  void group::accept(space &s) {
    int i,ilen=end+1;
    for (i=beg; i<ilen; i++) {
      s.p[i].x = s.trial[i].x;
      s.p[i].y = s.trial[i].y;
      s.p[i].z = s.trial[i].z;
    }
  }

  /*!
   * \param par Container class
   * \param c ...by adding this vector to all particles
   */
  void group::translate(container &con, const point &p) {
    for (int i=beg; i<=end; i++) {
      con.trial[i] = con.p[i] + p;
      con.boundary( con.trial[i] );
    }
  }
  
  /* -----------------------------------*
   *             MOLECULAR
   * -----------------------------------*/
   
  std::ostream& molecular::write(std::ostream &o) const {
    group::write(o);
    o << " " << cm;
    return o;
  }
  
  molecular& molecular::operator<<(std::istream &in) {
    group::operator<<(in);
    cm.operator<<(in);
    return *this;
  }
  
  void molecular::translate(container &con, const point &p) {
    group::translate(con, p);
    cm_trial=cm+p;
    con.boundary(cm_trial);
  }
  
  void molecular::accept(space &s) {
    group::accept(s);
    cm=cm_trial;
  }
  
}//namespace
