#include "faunus/hardsphere.h"

namespace Faunus {

  bool hardsphere::overlap(const vector<particle> &p, int j) {
    int ps=p.size();
    for (int i=0; i<j; i++)
      if (p[i].overlap(p[j])==true) return true;
    for (int i=j+1; i<ps; i++)
      if (p[i].overlap(p[j])==true) return true;
    return false;
  }

  bool hardsphere::overlap(const vector<particle> &p, const particle &a) {
    int i,n=p.size();
    for (i=0; i<n; i++)
      if (p[i].overlap(a)==true) return true;
    return false;
  }

  bool hardsphere::overlap(const vector<particle> &p, const group &g) {
    short int n=g.beg, psize=p.size();
    for (short int i=0; i<n; i++)
      if (overlap(p, g, i)==true) return true;
    for (short int i=n+1; i<psize; i++) 
      if (overlap(p, g, i)==true) return true;
    return false;
  }

  bool hardsphere::overlap(const vector<particle> &p, const group &g, int j) {
    if (g.beg==-1) return false;
    if (g.cm.radius>0)
      if ( abs( sqrt(g.cm.sqdist(p[j])) ) > g.cm.radius+p[j].radius)
        return false;
    int len=g.end+1;

    if (g.find(j)==false) {          //check if j is part of g
      for (int i=g.beg; i<len; i++)
        if (p[i].overlap(p[j])==true)
          return true;
    } else {                              //it is! avoid self-overlap...
      for (int i=g.beg; i<j; i++)
        if (p[i].overlap(p[j])==true) return true;
      for (int i=j+1; i<len; i++)
        if (p[i].overlap(p[j])==true) return true;
    }
    return false;
  }

  bool hardsphere::overlap(const vector<particle> &p, const group &g, const particle &a) {
    int len=g.end+1;
    if (g.beg!=-1)
      for (int i=g.beg; i<len; i++)
        if (p[i].overlap(a)==true)
          return true;
    return false;
  }

  bool hardsphere::overlap(const vector<particle> &p, const group &g1, const group &g2) {
    if (g1.beg==-1 || g2.beg==-1) return false;
    if (g1.cm.radius!=0 && g2.cm.radius!=0)
      if (abs( sqrt(g1.cm.sqdist(g2.cm)) ) > g1.cm.radius+g2.cm.radius)
        return false;
    int ilen=g1.end+1, jlen=g2.end+1;
    for (int i=g1.beg; i<ilen; i++) {
      for (int j=g2.beg; j<jlen; j++) {
        if ( p[i].overlap(p[j])==true )
          return true;
      }
    }
    return false;
  }
  bool hardsphere::overlap(const vector<particle> &p, const group &g1, const group &g2, double &s) {
    if (g1.beg==-1 || g2.beg==-1) return false;
    if (g1.cm.radius!=0 && g2.cm.radius!=0)
      if (abs( sqrt(g1.cm.sqdist(g2.cm)) ) > g1.cm.radius+g2.cm.radius)
        return false;
    int ilen=g1.end+1, jlen=g2.end+1;
    for (int i=g1.beg; i<ilen; i++) {
      for (int j=g2.beg; j<jlen; j++) {
        if ( p[i].overlap(p[j], s)==true )
          return true;
      }
    }
    return false;
  }
}//namespace
