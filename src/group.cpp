#include <faunus/geometry.h>
#include <faunus/point.h>
#include <faunus/species.h>
#include <faunus/physconst.h>
#include <faunus/group.h>
#include <faunus/average.h>
#include <faunus/space.h>
#include <faunus/inputfile.h>
#include <faunus/textio.h>

namespace Faunus {

  Group::Group(int front, int back) : Range(front,back-front+1) {
    if (front<0 || back<0)
      resize(0);
  }

  Group::~Group() {}

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

  string Group::_info() {
    std::ostringstream o;
    return o.str();
  }

  string Group::info() {
    using namespace textio;
    char w=15;
    std::ostringstream o;
    o << header("Group: " + name)
      << pad(SUB,w,"Size") << size() << endl;
    if (!empty())
      o << pad(SUB,w,"Range") << "[" << front() << "-" << back() << "]" << endl;
    o << pad(SUB,w,"Mass center") << cm.transpose() << endl;
    return o.str() + _info();
  }

  void Group::scale(Space &s, double newvol) {
    if (!empty()) {
      cm_trial=cm;
      cm_trial.scale(*s.geo, newvol);
      for (auto i : *this)
        s.trial[i].scale(*s.geo, newvol);
    }
  }

  int Group::random() const {
    if (!empty()) {
      int i = front() + slp_global.rand() % size();
      assert(find(i) && "Generated random element out of range!");
      return i;
    }
    return -1;
  }

  bool Group::isAtomic() const { return false; }

  bool Group::isMolecular() const { return false; }

  int Group::numMolecules() const {
    return size();
  }

  GroupAtomic::GroupAtomic(int front, int back) : Group(front, back) {
  }

  GroupAtomic::GroupAtomic(Space &spc, InputMap &in) {
    add(spc,in);
  }

  bool GroupAtomic::isAtomic() const { return true; }

  GroupMolecular::GroupMolecular(int front, int back) : Group(front, back) {
  }

  bool GroupMolecular::isMolecular() const { return true; }

  int GroupMolecular::numMolecules() const { return 1; }

  void GroupMolecular::scale(Space &s, double newvol) {
    assert( s.geo->dist(cm, massCenter(s))<1e-6);
    assert( s.geo->dist(cm, cm_trial)<1e-7);

    Point newcm=cm;
    newcm.scale(*s.geo, newvol);
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

  string GroupMolecular::_info() {
    using namespace textio;
    std::ostringstream o;
    return o.str();
  }

  GroupArray::GroupArray(int AtomsPerMolecule) : N(AtomsPerMolecule) {
    assert(N>0 && "Molecular size must be larger than zero");
  }

  int GroupArray::randomMol() const {
    int i=(random()-front())/N;
    assert(find(i) && "Out of range!");
    return i;
  }

  int GroupArray::numMolecules() const {
    assert( (size()%N)==0 );
    return size()/N;
  }

  /**
   * @warning You must manually update the mass center of the returned group
   */
  GroupMolecular& GroupArray::operator[](int i) {
    sel.setfront( front()+i*N );
    sel.setback( sel.front()+N-1 );
    assert( sel.size()==N );
    assert( find( sel.front() ) );
    assert( find( sel.back()  ) );
    return sel;
  }

  string GroupArray::_info() {
    using namespace textio;
    char w=15;
    std::ostringstream o;
    o << pad(SUB,w,"Mol size") << N << endl
      << pad(SUB,w,"Molecules") << numMolecules() << endl;
    return o.str();
  }

  void GroupArray::add(const GroupMolecular &g) {
    if ((g.size()%N)==0) {
      if (empty())
        setrange(g.front(), g.back());
      else if (g.front()==back()+1)
        setback(g.back());
    }
    assert( (size()%N)==0 && "GroupArray not a multiple of N");
  }

  void GroupArray::scale(Space &s, double newvol) {
    for (int i=0; i<numMolecules(); i++) {
      operator[](i);
      sel.setMassCenter(s);
      sel.scale(s,newvol);
    }
  }

}//namespace
