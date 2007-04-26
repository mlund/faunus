#include "group.h"

group::group(int i) {
  set(-1,-1);     //is empty
  chain=false;    //is not a chain
  graftpoint=-1;  //is not grafted
  radius=0;       //is sizeless
  name="(noname)";//is nameless
  vdw=false;      //omit vdw interactions
}

void group::set(short int first, short int last) {
  beg=first;
  end=last;
}

short int group::size() { return (beg==-1) ? 0 : end-beg+1; }
void group::operator++ (int) { end++; }
void group::operator-- (int) { end--;}

//append another group.
//they 
void group::operator+=(group g) {
  if (g.beg==-1 && g.end==-1)
    return;
  if (beg==-1 && end==-1)
    beg=g.beg;
  end=g.end;
}

//merge two groups into one. Their ranges must be continous...!
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

bool group::find(unsigned int i)
{ return (i>=beg && i<=end) ? true : false; }

