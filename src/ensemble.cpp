#include <faunus/ensemble.h>
#include "faunus/species.h"
#include "faunus/group.h"
#include "faunus/point.h"
#include "faunus/container.h"

namespace Faunus {

  //
  // class ensemble
  //

  ensemble::ensemble() {}

  bool ensemble::metropolis(double du)
  {
    if (du > 0)
      if ( slp.random_one()>exp(-du) )
        return false;
    return true;
  }

  //
  // class canonical
  //

  canonical::canonical() {}

  //
  // class isobarical
  //

  isobarical::isobarical() {}

  //
  // class grandcanonical
  //
  string grandcanonical::print(){
    inputfile in;
    int conpmax=0;
    for (int i=0; i<gp.size(); i++) {
      if (gp[i]->end >conpmax)
        conpmax=gp[i]->end;
      in.add(gp[i]->name+"beg", gp[i]->beg);
      in.add(gp[i]->name+"end", gp[i]->end);
    }
    in.add("consize", conpmax+1);
    return in.print();
  }
  // Loads correct size of groups and scales the particle
  // vector to continue a simulation. Call this before 
  // ioaam::load
  bool grandcanonical::load(container &con, string name) {
    inputfile in;
    if(in.load(name)==false)
      return false;
    for (int i=0; i<gp.size(); i++) {
      gp[i]->beg=in.getint(gp[i]->name+"beg",-1);
      gp[i]->end=in.getint(gp[i]->name+"end",-1);
    }
    con.p.clear();
    con.trial.clear();
    con.p.resize(in.getint("consize",-1));
    con.p.resize(in.getint("consize",-1));
    return true;
  }

  int grandcanonical::size() { return 0; }

  short grandcanonical::findgroup(unsigned int i) {
    for (unsigned int n=0; n<gp.size(); n++)
      if (gp[n]->find(i)==true) return n;
    return -1;
  }

  void grandcanonical::addGroupPtr(group &group) {
    gp.push_back(&group);
  }

  short grandcanonical::findgroup(string name) {
    for (int k=0; k<gp.size(); k++) {
      if (gp[k]->name==name) 
        return k;
    }
    return -1;
  }

  //\param p - particle vector to insert into.
  //\param i - We will insert at this position - the rest will be moved forward.
  //\param a - Particle to insert
  //\param n - index to group with group.name==atom[a.id].name
  // NOTE: These functions assume that any groups for ionic specis with chemical
  //       potential are NOT overlaping. Any other groups will be inplicitly handled
  //       by grand canonical rutines, overlaping or not. The functions not fully general.
  //       Will not handle groups with endpoints in between endpoints of saltgroups properly.
  bool grandcanonical::insert(container &con, particle &a) {
    short n=findgroup(atom[a.id].name);
    short i,j;
    i=gp[n]->beg, j=gp[n]->end;
    if (n>-1) {
      con.insert( a, gp[n]->end+1);
      for (unsigned short k=0; k<gp.size(); k++) {
        if (i<gp[k]->beg) 
          gp[k]->beg ++;
        if (j<=gp[k]->end)
          gp[k]->end ++;
      }  
      return true;
    }
    return false;
  }

  bool grandcanonical::erase(container &con, unsigned int &i) {
    short n=findgroup(atom[con.p[i].id].name);
    short l,j;
    l=gp[n]->beg, j=gp[n]->end;
    if (n>-1) {
      for (unsigned short k=0; k<gp.size(); k++) {
        if (l<gp[k]->beg) 
          gp[k]->beg --;
        if (j<=gp[k]->end)
          gp[k]->end --;
        }
      con.remove( i );
      return true;
    }
    return false;
  }

  bool grandcanonical::meMbeR(unsigned &i, unsigned int &j) {
    if (gp[i]->end>j && gp[i]->beg<=j)
      return true;
    else
      return false;
  }

  string grandcanonical::info() {
    std::ostringstream o;
    o << "Grand Canonical Master Class:\n"
      << "  Number of *groups     = " << gp.size() << endl
      << "  Number of particles  = " << size() << endl;
    return o.str();
  }

} //namespace
