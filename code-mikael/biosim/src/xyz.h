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
    void writetrj( vector<particle> &p, peptide &pep ) {
      f << p.size() << endl;
      if (virgin==true) f << endl; // VMD hack
      for (int i=0; i<p.size(); i++)
        f << pep.d[p[i].id].name << " "
          << p[i].x << " " << p[i].y << " " << p[i].z << endl;
      f << endl;
      virgin=false;
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


#endif
