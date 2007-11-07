#include "group.h"

group::group(int i) {
  set(-1,-1);     //is empty
  chain=false;    //is not a chain
  charge=0;       //is neutral
  graftpoint=-1;  //is not grafted
  cm=-1;          //is weightless
  radius=0;       //is sizeless
  name="(noname)";//is nameless
  dipv=-1;        //has no dipole moment
  vdw=false;      //omit vdw interactions
}

void group::set(short int first, short int last) {
  beg=first;
  end=last;
}

short int group::size() {
  if (beg==-1)
    return 0;
  return end-beg+1;
}
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
    cout << "# Warning: Added groups not continous!\n";
  return o;
}

short int group::random() { return beg + rand() % size(); }

ostream &operator<<( ostream &out, group &g) {
  out << "# " << g.name << endl
      << "#   Particles     = " << g.size() << " [" << g.beg << "," << g.end << "]\n"
      << "#   Charge        = " << g.charge << endl;
  if (g.graftpoint>=0)
    out<<"#   Graftpoint    = " << g.graftpoint << endl;
  if (g.radius>0)
    out<<"#   Radius        = "<< g.radius << endl;
  if (g.dipv>-1)
    out<<"#   Dipole vector = "<< g.dipv << endl;
  if (g.cm>-1)
    out<<"#   CM vector     = "<< g.cm << endl;

  return out;
}

bool group::isingroup(int i) {
  if (i>=beg && i<=end)
    return true;
  return false;
}
