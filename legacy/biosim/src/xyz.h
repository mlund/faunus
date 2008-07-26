#ifndef _XYZ_H
#define _XYZ_H

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include "point.h"
#include "peptide.h"

class xyz {
  private:
    unsigned short fid;
    ofstream f;
    string base;
    bool virgin;
  public:
    xyz(string file="trj") {
      fid=0;
      base=file;
      newfile();
    }
    void header( unsigned short n ) {
      f << n << endl;
      if (virgin==true) f << endl;
    }
    void footer() { f << endl; }
    void add( vector<particle> &p, peptide &pep, unsigned short i) {
      f << pep.d[p[i].id].name.substr(0,3) << " "
                  << p[i].x << " " << p[i].y << " " << p[i].z << endl;
      virgin=false;
    }
    void add( vector<particle> &p, peptide &pep, group &g) {
      for (unsigned short i=g.beg; i<=g.end-2; i++) add(p, pep, i); //skip cm and dipole!
    }
    void newfile() {
      if (f) close();
      ostringstream suffix;
      suffix << "." << ++fid << ".xyz";
      string s = base + suffix.str();
      f.open( s.c_str() );
      f.precision(3);
      virgin=true;
    }
    void close() { f.close(); }
};

// PQR file output stolen from faunus/classes.
class pqr {
  private:
    int nres, natom;
    string name;
    char buf[100];
    ofstream f;
  public:
    void add(vector<particle> &p, peptide &pep, group &g) {
      for (unsigned short i=g.beg; i<=g.end-2; i++) add(p,pep,i); //skip cm and dipole!
    }
    void add(vector<particle> &p, peptide &pep, unsigned short i) {
      name=pep.d[p[i].id].name.substr(0,3);
      sprintf(buf, "ATOM  %5d %-4s %s %5d    %8.3f %8.3f %8.3f %.3f %.3f\n",
          natom++,name.c_str(),name.c_str(),nres,p[i].x,p[i].y,p[i].z,p[i].charge,p[i].radius);
      f << buf;
      if (name.find("CTR")!=string::npos) nres++;
    }
    pqr(string filename) {
      nres=natom=1;
      f.open(filename.c_str());
      f.precision(3);
    }
    void close() { f.close(); }
};

#endif
