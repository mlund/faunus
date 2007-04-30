#include "group.h"

//---------------- GROUP -----------------------
group::group(int i) {
  set(-1,-1);     //is empty
  name="(noname)";//is nameless
}

void group::set(short int first, short int last) {
  beg=first;
  end=last;
}

short int group::size() { return (beg==-1) ? 0 : end-beg+1; }

void group::operator+=(group g) {
  if (g.beg==-1 && g.end==-1)
    return;
  if (beg==-1 && end==-1)
    beg=g.beg;
  end=g.end;
}

group group::operator+(group g) {
  group o;
  if (this->beg<g.beg) {
    o.beg=this->beg;
    o.end=g.end;
  } else {
    o.beg=g.beg;
    o.end=this->end;
  };
  o.name = this->name + " + " + g.name;
  if (o.size()!=this->size()+g.size())
    cout << "# Warning: Added groups are not continous!\n";
  return o;
}

short int group::random() { return beg + rand() % size(); }

ostream &operator<<( ostream &out, group &g) {
  stringstream s;
  s << "[" << g.beg << "," << g.end << "]";
  //s.str(""); // clear stringstream
  out << "# "
    << left
    << setw(10) << g.name
    << right
    << setw(10) << s.str()
    << setw(5)  << g.Q.avg()
    << setw(22) << g.cm
    << setw(22) << g.mu
    << endl;
  return out;
}

bool group::find(unsigned int i) { return (i>=beg && i<=end) ? true : false; }

//--------------- MACRO MOLECULE ---------------
macromolecule::macromolecule() {
  radius=0;       //is sizeless
}


//--------------- CHAIN -----------------
chain::chain() {
  graftpoint=-1;  //is not grafted
}
