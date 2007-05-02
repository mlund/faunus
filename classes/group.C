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

/*! \brief Calculates total charge
 *  \return Charge number. The charge is also stored in
 *          cm.charge.
 */
double group::charge(vector<particle> &p) {
  double z=0;
  for (short i=beg; i<=end; i++)
    z+=p[i].charge;
  cm.charge=z;
  return z;
}

void group::operator+=(group g) {
  if (g.beg==-1 && g.end==-1)   // if added group is empty
    return;
  if (beg==-1 && end==-1)       // if this is empty
    beg=g.beg;
  end=g.end;
}

group group::operator+(group g) {
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
    << setw(22) << g.cm
    << endl;
  return out;
}

bool group::find(unsigned int i) { return (i>=beg && i<=end) ? true : false; }

/*!
 * Recalcute mass center of the particles
 * in a group. The result is stored in the
 * group CM placeholder as well as returned
 * by the function.
 */
point group::masscenter(vector<particle> &p) {
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

//--------------- MACRO MOLECULE ---------------
macromolecule::macromolecule() {}
void macromolecule::operator=(group g) {
  beg=g.beg;
  end=g.end;
  name=g.name;
  cm=g.cm;
  cm_trial=g.cm_trial;
}
/*!
 * Recalulates the dipole moment of a
 * group and stores the result in the
 * group dipole moment placeholder.
 * Origo is (0,0,0).
 */
double macromolecule::dipole(vector<particle> &p)
{
  mu.clear();
  for (int i=beg; i<=end; i++) {
    mu.x += (p[i].x-cm.x) * p[i].charge;
    mu.y += (p[i].y-cm.y) * p[i].charge;
    mu.z += (p[i].z-cm.z) * p[i].charge;
  }
  return mu.len();
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


//--------------- CHAIN -----------------
chain::chain() {
  graftpoint=-1;  //is not grafted
}

double chain::monomerenergy(vector<particle> &p, short i) {
  double u=0;
  //the first ?
  if (i==beg) {
    u+=quadratic( p[i], p[i+1] );
    if (graftpoint>-1)
      u+=quadratic( p[i], p[graftpoint] );
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
  if (graftpoint>-1)
    u+=quadratic( p[beg], p[graftpoint] );
  return u;
}
